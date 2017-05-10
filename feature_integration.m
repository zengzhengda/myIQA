%% feature integration when some features are changed
clear ;
clc;
close all;

database =1; % 1:CID2013, 2:LIVE, 3:TID2013
if(database==1)
    load('my_mat_cid2013_salient.mat');
    salient_fea=my_mat;
    load('my_mat_cid2013.mat');
    
    start=2;
    num1=75;
    num_fea=(size(my_mat,2)-2)/2;
    my_mat=[my_mat(:,1) my_mat(:,start:start+num1-1) salient_fea(:,2:4) my_mat(:,start+num_fea : start+num_fea+num1-1) salient_fea(:,5:7) my_mat(:,end)];
    save my_mat_cid2013_20170507.mat my_mat;   
elseif(database==2)
    load('my_mat_live_salient.mat');
    salient_fea=my_mat;
    load('my_mat_live.mat');
    
    start=2;
    num1=75;
    num_fea=(size(my_mat,2)-2)/2;
    my_mat=[my_mat(:,1) my_mat(:,start:start+num1-1) salient_fea(:,2:4) my_mat(:,start+num_fea : start+num_fea+num1-1) salient_fea(:,5:7) my_mat(:,end)];
    save my_mat_live_20170506.mat my_mat; 
else
    load('my_mat_tid2013_salient.mat');
    salient_fea=my_mat;
    load('my_mat_tid2013.mat');
    
    start=2;
    num1=75;
    num_fea=(size(my_mat,2)-2)/2;
    my_mat=[my_mat(:,1) my_mat(:,start:start+num1-1) salient_fea(:,2:4) my_mat(:,start+num_fea : start+num_fea+num1-1) salient_fea(:,5:7) my_mat(:,end)];
    save my_mat_tid2013_20170506.mat my_mat; 
end