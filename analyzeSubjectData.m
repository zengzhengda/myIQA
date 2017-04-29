% �����������ݼ�
clear;
clc;
path_excel='..\Datasets\CID2013\CID2013 data - version 12112014.xlsx';
[number,txt,raw]=xlsread(path_excel,'CID2013 MOS');
realigned_mos=number(:,3);
sharpness=number(:,4);
graininess=number(:,5);
lightness=number(:,6);
color_saturation=number(:,9);
feature=[sharpness graininess lightness];
subject_data=[feature realigned_mos];

%% ��׼������
mu_my=mean(subject_data(:,1:end));
sigma2_my=mean((subject_data(:,1:end)-repmat(mu_my,474,1)).^2);
subject_data=(subject_data(:,1:end)-repmat(mu_my,474,1))./repmat((sigma2_my.^0.5),474,1);
%% ѵ�����ݺͲ��������и�
cnt_groups=6;
cnt_train_groups=4;
[train_data,test_data]=crossValidation(subject_data,cnt_groups,cnt_train_groups);
    
%% svrģ��ѵ��
% ѡ�����Ų���
[best_C,best_gamma]=SVR_choosing_paremeter(train_data(:,1:end-1),train_data(:,end));
model=svmtrain(train_data(:,end),train_data(:,1:end-1),sprintf('-s %f -t %f -c %f -g %f', 3, 2, best_C, best_gamma)); % add -v ������֤�������ģ��ֵ
%% ѵ������֤���
[train_quality_norm,train_mse,train_decision]=svmpredict(train_data(:,end),train_data(:,1:end-1),model);
train_quality=train_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
train_mos=train_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%% svr Ԥ��
[test_quality_norm,test_mse,test_decision]=svmpredict(test_data(:,end),test_data(:,1:end-1),model);
test_quality=test_quality_norm.*(sigma2_my(:,end)^0.5)+mu_my(:,end);
test_mos=test_data(:,end).*(sigma2_my(:,end)^0.5)+mu_my(:,end);
%% ��������
isShow=1;
[CC,SROCC,RMSE]=performance_eval(test_quality,test_mos,isShow);