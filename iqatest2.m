clear ;
clc;
close all;
warning('off'); % ����ʾwarning
                    
% ѡ�����ݿ�
dataset_names={'CID2013','LIVE','TID2013'};
dataset_name=dataset_names{1};% ѡ�����ݿ�
% metric choice
metric='Robust_CIQA';  % Robust_CIQA,Robust_CIQA_NoNorm,IDEAL,BRISQUE,BLIINDS2

% load best_svr_param_tid2013;

%% ������ȡ   
% fetchFeatureAll(dataset_name);
%% ��������
% load(['my_mat_' lower(dataset_name)]);
load('my_mat_cid2013');
% load('my_mat_tid2013_20170517');
% ȥ���ο�ͼ��
%  load('..\Datasets\LIVE\dmos_realigned.mat');
%  orgs2=orgs(1:end-174);
% inx=find(orgs2==1);
%  my_mat(inx,:)=[];
%% ����ѡ��
if(strcmp(metric,'Robust_CIQA'))
    num_fea=(size(my_mat,2)-2)/2;
    num_base=6;
    num_grad=6;
    num_nss=25;
    num_sharp=1;
    num_salient=2;
    % ��������
    start=14; % 2��һ��ģ�14��³����
    base_fea=[my_mat(:,start:start+num_base-1) my_mat(:,start+num_fea : start+num_base+num_fea-1)];
    %�ݶ�����
    start=20; % 8��һ��ģ�20��³����
    grad_fea=[my_mat(:,start:start+num_grad-1) my_mat(:,start+num_fea : start+num_grad+num_fea-1)];
    % nss���� 
    start=51; % 26��һ��ģ�51��³����
    nss_fea=[my_mat(:,start:start+num_nss-1) my_mat(:,start+num_fea : start+num_nss+num_fea-1)];
    % sharp����
    start=76;
    sharp_fea=[my_mat(:,start:start+num_sharp-1) my_mat(:,start+num_fea : start+num_sharp+num_fea-1)];
    % salient ����
    start=77;
    salient_fea=[my_mat(:,start:start+num_salient-1) my_mat(:,start+num_fea : start+num_salient+num_fea-1)];
 
    my_mat=[my_mat(:,1) base_fea grad_fea nss_fea sharp_fea salient_fea my_mat(:,end)]; 
end
%% ȥ�������ʵ�����
my_mat(404,:)=[];
num_data=size(my_mat,1);

%% ������֤
cnt_cross=1000;% ������֤����

CC_my=zeros(1,cnt_cross);
SROCC_my=zeros(1,cnt_cross);
RMSE_my=zeros(1,cnt_cross);

CC_ideal=zeros(1,cnt_cross);
SROCC_ideal=zeros(1,cnt_cross);
RMSE_ideal=zeros(1,cnt_cross);

CC_brisque=zeros(1,cnt_cross);
SROCC_brisque=zeros(1,cnt_cross);
RMSE_brisque=zeros(1,cnt_cross);


isBestParam=0;% ��ʶ�����Ƿ�Ϊsvm���Ų���
isShow=1;  % �Ƿ���ʾ����Ա�ͼ

rate_train=0.8; % ѵ��������
rate_test=0.2;  % ���Լ�����

for ii=1:cnt_cross
    [train_inds,test_inds]=randomSample(num_data,rate_train,rate_test);
    switch metric
        case 'Robust_CIQA'         
            [train_data,test_data]=dataSplit(my_mat,train_inds,test_inds);
            [CC_my(ii),SROCC_my(ii),RMSE_my(ii)]=robustSVR(train_data,test_data,dataset_name,isBestParam,isShow);          
            isShow=0;
            isBestParam=1;
%             start=2; % �������
%             end2=1; % �����յ�
% 
%             %% ��׼������ ���ԸĽ�
%             mu_my=mean(my_mat(:,2:end));
%             sigma2_my=mean((my_mat(:,2:end)-repmat(mu_my,num_data,1)).^2);
%             my_mat_norm=[my_mat(:,1),(my_mat(:,2:end)-repmat(mu_my,num_data,1))./repmat((sigma2_my.^0.5),num_data,1)];
% 
%             [train_data,test_data]=dataSplit(my_mat_norm,train_inds,test_inds);
%             %% svrģ��ѵ��
%             % ��������������ţ���ѡ�����Ų���
%             if(isBestParam==0)
%                 [best_C,best_gamma]=SVR_choosing_paremeter(train_data(:,start:end-end2),train_data(:,end),dataset_name);
%                 isBestParam=1;
%             end
%             model=svmtrain(train_data(:,end),train_data(:,start:end-end2),sprintf('-s %f -t %f -c %f -g %f', 3, 2, best_C, best_gamma)); % add -v ������֤�������ģ��ֵ
%             %% ѵ������֤���
%             [train_quality_norm,train_mse,train_decision]=svmpredict(train_data(:,end),train_data(:,start:end-end2),model);
%             train_quality=train_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%             train_mos=train_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%             
%             %% svr Ԥ��
%             [test_quality_norm,test_mse,test_decision]=svmpredict(test_data(:,end),test_data(:,start:end-end2),model);
%             test_quality=test_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%             test_mos=test_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%             %% ��������
%             [CC_my(ii),SROCC_my(ii),RMSE_my(ii)]=performance_eval(test_quality,test_mos,isShow);
%             isShow=0;
        case 'Robust_CIQA_NoNorm'  % û�б�׼������
            start=2; % �������
            end2=1; % �����յ�
            [train_inds,test_inds]=randomSample(num_data,rate_train,rate_test);
            [train_data,test_data]=dataSplit(my_mat,train_inds,test_inds);
            if(isBestParam==0)
                [best_C,best_gamma]=SVR_choosing_paremeter(train_data(:,start:end-end2),train_data(:,end),dataset_name);
                isBestParam=1;
            end
            model=svmtrain(train_data(:,end),train_data(:,start:end-end2),sprintf('-s %f -t %f -c %f -g %f', 3, 2, best_C, best_gamma)); % add -v ������֤�������ģ��ֵ
            %% svr Ԥ��
            [test_quality,test_mse,test_decision]=svmpredict(test_data(:,end),test_data(:,start:end-end2),model);
            %% ��������
            test_mos=test_data(:,end);
            [CC_my(ii),SROCC_my(ii),RMSE_my(ii)]=performance_eval(test_quality,test_mos,isShow);
            isShow=0;
        case 'BRISQUE'
            
    end
   
end
mu_CC=mean(CC_my);
mu_SROCC=mean(SROCC_my);
mu_RMSE=mean(RMSE_my);
med_CC=median(CC_my);
med_SROCC=median(SROCC_my);
med_RMSE=median(RMSE_my);
result_mean=['|' num2str(mu_CC) '|' num2str(mu_SROCC) '|'  num2str(mu_RMSE) '|']
result_med=['|' num2str(med_CC) '|' num2str(med_SROCC) '|'  num2str(med_RMSE) '|']
    