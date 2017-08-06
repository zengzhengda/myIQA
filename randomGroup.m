function [data_cell]=randomGroup(data,groups)
cnt_data=size(data,1);
cnt_per_group=floor(cnt_data/groups); % the number of each group
i=0;

data_cell={};
for i =1:groups-1
	rng('shuffle');
	rand_indexs=randperm(cnt_data-cnt_per_group*(i-1),cnt_per_group);
	data_cell{i}=data(rand_indexs,:);
    data(rand_indexs,:)=[];
end
data_cell{i+1}=data;