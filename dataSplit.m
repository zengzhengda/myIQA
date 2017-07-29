function [train_data,test_data]=dataSplit(data,train_inds,test_inds)
train_data=data(train_inds,:);
test_data=data(test_inds,:);
end

