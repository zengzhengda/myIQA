function [CC_my,SROCC_my,RMSE_my]=robustSVR(train_data,test_data,dataset_name,isBestParam,isShow)
start=2; % 特征起点
end2=1; % 特征终点

%% 标准化处理 可以改进
num_train=size(train_data,1);
num_test=size(test_data,1);

mu_my=mean(train_data(:,2:end));
sigma2_my=mean((train_data(:,2:end)-repmat(mu_my,num_train,1)).^2);
train_data=[train_data(:,1),(train_data(:,2:end)-repmat(mu_my,num_train,1))./repmat((sigma2_my.^0.5),num_train,1)];
test_data=[test_data(:,1),(test_data(:,2:end)-repmat(mu_my,num_test,1))./repmat((sigma2_my.^0.5),num_test,1)];

%% svr模型训练
% 如果参数不是最优，则选择最优参数
if(isBestParam==0)
    [best_C,best_gamma]=SVR_choosing_paremeter(train_data(:,start:end-end2),train_data(:,end),dataset_name);
else
    load best_svr_param_cid.mat
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
[CC_my,SROCC_my,RMSE_my]=performance_eval(test_quality,test_mos,isShow);
