function [train_inds,test_inds]=randomSample(num_data,rate_train,rate_test)
num_train=floor(num_data*rate_train);
rng('shuffle')
data_inds=[1:num_data];
train_inds=randperm(num_data,num_train);
data_inds(train_inds)=[];
test_inds=data_inds;
end