function []=trainOfCluter(img_ind,train_quality)
%% ����train_quality ��ͼ�񼯷���
cnt_train=length(train_quality);
level_vec=zeros(cnt_train,1); %  �ֳ�10��
for i=1:cnt_train
    level_vec(i)=ceil(train_quality(i)/10);
end
img_info=[img_ind,level_vec];

end
