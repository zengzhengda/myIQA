clear all;
clc;
warning('off'); % 不显示warning
%% 特征提取                       
%CID2013
% fetchFeatureAll();
%% 性能评价
load my_mat;
%% 标准化处理
mu_my=mean(my_mat);
sigma2_my=mean((my_mat-repmat(mu_my,474,1)).^2);
my_mat=(my_mat-repmat(mu_my,474,1))./repmat((sigma2_my.^0.5),474,1);

%% 交叉验证
cnt_cross=100;% 交叉验证次数
CC=zeros(1,cnt_cross);
SROCC=zeros(1,cnt_cross);
RMSE=zeros(1,cnt_cross);
for ii=1:cnt_cross
    cnt_groups=5;
    cnt_train_groups=4;
    [train_data,test_data]=crossValidation(my_mat,cnt_groups,cnt_train_groups);
    %%
    model=svmtrain(train_data(:,end),train_data(:,1:end-1),'-s 3 -t 2 -c 2 -g 0.2 -p 0.01'); % add -v 交叉验证，不输出模型值
    [test_quality_norm,mse,decision]=svmpredict(test_data(:,end),test_data(:,1:end-1),model);
    test_quality=test_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
    test_mos=test_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
    %% 性能评估
    [CC(ii),SROCC(ii),RMSE(ii)]=performance_eval(test_quality,test_mos);
end
mu_CC=mean(CC)
mu_SROCC=mean(SROCC)
mu_RMSE=mean(RMSE)
        
    

