function [inds_cell]=getIndexsOfRandomGroup(cnt_data,groups)

cnt_per_group=floor(cnt_data/groups); % the number of each group
i=0;
inds_cell={};
inds_data=[1:cnt_data];
for i =1:groups-1
	rng('shuffle');  
	rand_indexs=randperm(cnt_data-cnt_per_group*(i-1),cnt_per_group);
	inds_cell{i}=inds_data(rand_indexs);
	inds_data(rand_indexs)=[];
end
inds_cell{i+1}=inds_data;