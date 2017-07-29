clear ;
clc;
close all;
warning('off'); % 不显示warning
                    
% 选择数据库
dataset_names={'CID2013','LIVE','TID2013'};
dataset_name=dataset_names{1};% 选择数据库
% load best_svr_param_tid2013;

%% 特征提取   
% fetchFeatureAll(dataset_name);
%% 性能评价
% load(['my_mat_' lower(dataset_name)]);
load('my_mat_cid2013');
% load('my_mat_tid2013_20170517');
% 去除参考图像
%  load('..\Datasets\LIVE\dmos_realigned.mat');
%  orgs2=orgs(1:end-174);
% inx=find(orgs2==1);
%  my_mat(inx,:)=[];
%% 特征选择
isFeatureAnalyse=0; % 是否单项分析
if(isFeatureAnalyse ~=1)
    num_fea=(size(my_mat,2)-2)/2;
    num_base=6;
    num_grad=6;
    num_nss=25;
    num_sharp=1;
    num_salient=2;
    % 基本特征
    start=14; % 2是一般的，14是鲁棒的
    base_fea=[my_mat(:,start:start+num_base-1) my_mat(:,start+num_fea : start+num_base+num_fea-1)];
    %梯度特征
    start=20; % 8是一般的，20是鲁棒的
    grad_fea=[my_mat(:,start:start+num_grad-1) my_mat(:,start+num_fea : start+num_grad+num_fea-1)];
    % nss特征 
    start=51; % 26是一般的，51是鲁棒的
    nss_fea=[my_mat(:,start:start+num_nss-1) my_mat(:,start+num_fea : start+num_nss+num_fea-1)];
    % sharp特征
    start=76;
    sharp_fea=[my_mat(:,start:start+num_sharp-1) my_mat(:,start+num_fea : start+num_sharp+num_fea-1)];
    % salient 特征
    start=77;
    salient_fea=[my_mat(:,start:start+num_salient-1) my_mat(:,start+num_fea : start+num_salient+num_fea-1)];
 
    my_mat=[my_mat(:,1) base_fea grad_fea nss_fea sharp_fea salient_fea my_mat(:,end)]; 
end
%% 去除不合适的数据
my_mat(404,:)=[];
cnt_data=size(my_mat,1);
%% 标准化处理 可以改进
mu_my=mean(my_mat(:,2:end));
sigma2_my=mean((my_mat(:,2:end)-repmat(mu_my,cnt_data,1)).^2);
my_mat=[my_mat(:,1),(my_mat(:,2:end)-repmat(mu_my,cnt_data,1))./repmat((sigma2_my.^0.5),cnt_data,1)];
%% 交叉验证
cnt_cross=1000;% 交叉验证次数
CC=zeros(1,cnt_cross);
SROCC=zeros(1,cnt_cross);
RMSE=zeros(1,cnt_cross);
isBestParam=0;% 标识符，是否为svm最优参数
isShow=1;  % 是否显示结果对比图

rate_train=0.8; % 训练集比例
rate_test=0.2;  % 测试集比例
num_data=size(my_mat,1);

% metric choice
metric='Robust_CIQA'
for ii=1:cnt_cross
    switch metric
        case 'Robust_CIQA'
            start=2; % 特征起点
            end2=1; % 特征终点
            [train_inds,test_inds]=randomSample(num_data,rate_train,rate_test);
            [train_data,test_data]=dataSplit(my_mat,train_inds,test_inds);
            %% svr模型训练
            % 如果参数不是最优，则选择最优参数
            if(isBestParam==0)
                [best_C,best_gamma]=SVR_choosing_paremeter(train_data(:,start:end-end2),train_data(:,end),dataset_name);
                isBestParam=1;
            end
            model=svmtrain(train_data(:,end),train_data(:,start:end-end2),sprintf('-s %f -t %f -c %f -g %f', 3, 2, best_C, best_gamma)); % add -v 交叉验证，不输出模型值
            %% 训练集验证结果
            [train_quality_norm,train_mse,train_decision]=svmpredict(train_data(:,end),train_data(:,start:end-end2),model);
            train_quality=train_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
            train_mos=train_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
            
            %% svr 预测
            [test_quality_norm,test_mse,test_decision]=svmpredict(test_data(:,end),test_data(:,start:end-end2),model);
            test_quality=test_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
            test_mos=test_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
            %% 性能评估
            [CC(ii),SROCC(ii),RMSE(ii)]=performance_eval(test_quality,test_mos,isShow);
            isShow=0;
        otherwise

    end
   
end
mu_CC=mean(CC);
mu_SROCC=mean(SROCC);
mu_RMSE=mean(RMSE);
med_CC=median(CC);
med_SROCC=median(SROCC);
med_RMSE=median(RMSE);
result_mean=['|' num2str(mu_CC) '|' num2str(mu_SROCC) '|'  num2str(mu_RMSE) '|']
result_med=['|' num2str(med_CC) '|' num2str(med_SROCC) '|'  num2str(med_RMSE) '|']
    