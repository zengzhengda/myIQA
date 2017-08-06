% 2017-07-31 test by cross validation
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
% load('my_mat_cid2013');
load('my_mat_live_20170517');
% ȥ���ο�ͼ��
%  load('..\Datasets\LIVE\dmos_realigned.mat');
%  orgs2=orgs(1:end-174);
% inx=find(orgs2==1);
%  my_mat(inx,:)=[];
%% ����ѡ��
if(strcmp(metric,'Robust_CIQA') || strcmp(metric,'Robust_CIQA_NoNorm'))
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
cnt_cross=6;% ������֤����

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
hasSampleIndex=0; % �Ƿ��������������

inds_cell={};

%% leave-one-out cross validation
best_C=0;
best_gamma=0;

for ii=1:cnt_cross
    if(hasSampleIndex ==0)
        inds_cell=getIndexsOfRandomGroup(num_data,cnt_cross);
        save indexs_randSample.mat inds_cell
    else
        load indexs_randSample.mat
    end
    inds_test=cell2mat(inds_cell(ii));
    inds_cell(ii)=[];
    inds_train=cell2mat(inds_cell);
    
    switch metric
        case 'Robust_CIQA'     
            train_data=my_mat(inds_train,:);
            test_data=my_mat(inds_test,:);      
            [CC_my(ii),SROCC_my(ii),RMSE_my(ii)]=robustSVR(train_data,test_data,dataset_name,isBestParam,isShow);          
            isShow=0;
            isBestParam=1;
        case 'Robust_CIQA_NoNorm'  % û�б�׼������
            start=2; % �������
            end2=1; % �����յ�
            train_data=my_mat(inds_train,:);
            test_data=my_mat(inds_test,:); 
            if(isBestParam==0)
                [best_C,best_gamma]=SVR_choosing_paremeter_no_norm(train_data(:,start:end-end2),train_data(:,end),dataset_name);
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
            
        case 'IDEAL'
%             start=2; % �������
%             end2=1; % �����յ�
%             load ideal_mat_cid2013.mat;
%             my_mat_ideal(404,:)=[];
%             train_data_ideal=my_mat_ideal(inds_train,:);
%             test_data_ideal=my_mat_ideal(inds_test,:);
%             train_DB='LIVE';
%             load(['Trained_SVR_module_on_' train_DB '.mat'],'model')
% 
%             [test_quality_ideal, ~, ~] = svmpredict(test_data_ideal(:,end),test_data_ideal(:,start:end-end2),model); 
%             test_mos_ideal=test_data_ideal(:,end);
%             [CC_ideal(ii),SROCC_ideal(ii),RMSE_ideal(ii)]=performance_eval(test_quality_ideal,test_mos_ideal,isShow);
%             isShow=0;

            start=2; % �������
            end2=1; % �����յ�
            load ideal_mat_cid2013.mat;
            my_mat_ideal(404,:)=[];
            train_data_ideal=my_mat_ideal(inds_train,:);
            test_data_ideal=my_mat_ideal(inds_test,:);
            if(isBestParam==0)
                [best_C,best_gamma]=SVR_choosing_paremeter_no_norm(train_data_ideal(:,start:end-end2),train_data_ideal(:,end),dataset_name);
                isBestParam=1;
            end
            model=svmtrain(train_data_ideal(:,end),train_data_ideal(:,start:end-end2),sprintf('-s %f -t %f -c %f -g %f', 3, 2, best_C, best_gamma)); % add -v ������֤�������ģ��ֵ
            %% svr Ԥ��
            [test_quality_ideal,test_mse_ideal,test_decision_ideal]=svmpredict(test_data_ideal(:,end),test_data_ideal(:,start:end-end2),model);
            %% ��������
            test_mos_ideal=test_data_ideal(:,end);
            [CC_ideal(ii),SROCC_ideal(ii),RMSE_ideal(ii)]=performance_eval(test_quality_ideal,test_mos_ideal,isShow);
            isShow=0;
    end
   
end
if(strcmp(metric,'Robust_CIQA'))
    showResults(SROCC_my,CC_my,RMSE_my,metric);
elseif(strcmp(metric,'IDEAL'))
    showResults(SROCC_ideal,CC_ideal,RMSE_ideal,metric);
elseif(strcmp(metric,'Robust_CIQA_NoNorm'))
    showResults(SROCC_my,CC_my,RMSE_my,metric);
end
