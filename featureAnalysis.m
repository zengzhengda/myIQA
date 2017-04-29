clc;
clear;
close all;

img={};
path{1}='IS_I_C03_D04_7.3705';
path{2}='IS_I_C03_D01_31.0333';
path{3}='IS_I_C03_D02_94.3333';

img_name{1}='IS\_I\_C03\_D04\_7.3705';
img_name{2}='IS\_I\_C03\_D01\_31.0333';
img_name{3}='IS\_I\_C03\_D02\_94.3333';

for i=1:length(path)
    img{i}=imread([path{i},'.jpg']);
end
cnt_img=length(img);
cnt_img_level=1;
for i=1:cnt_img
[fea{i} VSMap{i} grad_map1{i} structdis{i} M{i} N{i}]=fetchFeature2(img{i},cnt_img_level);
end
%% ������������ֱ��ͼͳ��
figure('name','������������ֱ��ͼͳ��');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    histogram(VSMap{i});
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});
%% ������������MSCN
figure('name','������������MSCN');
for i=1:cnt_img
mscn_VSMap{i}=MSCN(VSMap{i});
ind1=find(mscn_VSMap{i}>1);
ind2=find(mscn_VSMap{i}<-1);
mscn_VSMap{i}(ind1)=1;
mscn_VSMap{i}(ind2)=-1;
subplot(1,cnt_img,i);
h_mscn_vs{i}=histogram(mscn_VSMap{i});
legend(img_name{i});
hold on;
end


figure('name','����������MSCN���')
for i=1:cnt_img
    h_values=h_mscn_vs{i}.Values;
    h_values2=h_values(1:4:end);
    plot(h_values2,'*-');
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});
%% ԭͼ��ֱ��ͼ
% figure('name','ԭͼ��ֱ��ͼ');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(img{i}(:));
%     xlabel(img_name{i});
% end
%% M������ֱ��ͼ
% figure('name','M������ֱ��ͼ');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(M{i}(:));
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});



%% N������ֱ��ͼ
% figure('name','N������ֱ��ͼ');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(N{i}(:));
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
%% M������MSCN
% figure('name','M������MSCN');
% for i=1:cnt_img
% mscn_M{i}=MSCN(M{i});
% % subplot(1,cnt_img,i);
% h_mscn_M{i}=histogram(mscn_M{i});
% hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
% 
% figure('name','M������MSCN���')
% for i=1:cnt_img
%     h_values=h_mscn_M{i}.Values;
%     h_values=medfilt1(h_values,3);
%     plot(h_values,'-');
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
%% N������MSCN
% figure('name','N������MSCN');
% for i=1:cnt_img
% mscn_N{i}=MSCN(M{i});
% % subplot(1,cnt_img,i);
% h_mscn_N{i}=histogram(mscn_N{i});
% hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
% 
% figure('name','N������MSCN���')
% for i=1:cnt_img
%     h_values=h_mscn_N{i}.Values;
%     plot(h_values,'-');
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});

%% ��MSCN��ֱ��ͼ
% figure('name','��MSCN��ֱ��ͼ');
% for i=1:cnt_img
%     ind1=find(structdis{i}>2);
%     ind2=find(structdis{i}<-2);
%     structdis{i}(ind1)=2;
%     structdis{i}(ind2)=-2;
% %     subplot(1,cnt_img,i);
%     histogram(structdis{i});
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});