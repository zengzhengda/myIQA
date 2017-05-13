clc;
clear;
close all;

img={};
% path{1}='IS_I_C04_D06_54.3';
% path{2}='IS_IV_C04_D07_69.4';
% path{3}='IS_I_C04_D11_83.0';
% img_name{1}='MOS=54.3';
% img_name{2}='MOS=69.4';
% img_name{3}='MOS=83.0';

path{1}='IS_VI_C06_D01_29.8';
path{2}='IS_VI_C06_D05_64.8';
path{3}='IS_VI_C06_D10_87.9';
img_name{1}='MOS=29.8';
img_name{2}='MOS=64.8';
img_name{3}='MOS=87.9';

trunc_val=0.2e1;
tmp1=1e17;
thresh_y=5e5;

for i=1:length(path)
    img{i}=imread([path{i},'.jpg']);
end
cnt_img=length(img);
cnt_img_level=1;
for i=1:cnt_img
[structdis_L{i},structdis_M{i},structdis_N{i},structdis_Lmax{i},structdis_Lmin{i},feature]=analyseNSSFeature(img{i});
end
%% L的MDSCN特征的直方图统计
figure('name','L的MDSCN的直方图统计');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    ind1=find(structdis_L{i}>trunc_val);
    ind2=find(structdis_L{i}<-trunc_val);
    ind=[ind1;ind2];
    structdis_L{i}(ind)=[];
    h_L=histogram(structdis_L{i});
    h_values{i}=h_L.Values;
    h_edges{i}=h_L.BinEdges;
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});

%L的折线图
figure('name','L的MDSCN的折线图');
for i=1:cnt_img
    y=h_values{i};
    x=h_edges{i};
    ind=find(y~=0);
    y2=y(ind);
    x2=x(ind)/(tmp1);
    plot(x2,y2,'*-');
    hold on;
end
xlabel('MDSCN of Luminance')
ylabel('Number of Coffecients')
legend(img_name{1},img_name{2},img_name{3});

%% M的MDSCN特征的直方图统计
figure('name','M的MDSCN的直方图统计');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    structdis=structdis_M{i};
    ind1=find(structdis>trunc_val*0.01);
    ind2=find(structdis<-trunc_val*0.01);
    ind=[ind1;ind2];
    structdis(ind)=[];
    h=histogram(structdis);
    h_values{i}=h.Values;
    h_edges{i}=h.BinEdges;
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});

%% M的折线图
figure('name','M的MDSCN的折线图');
for i=1:cnt_img
    y=h_values{i};
    x=h_edges{i};
    [~,ind_max]=max(y);
    if(y(ind_max-1) < 5e4)
        y(ind_max-1)=[];
        x(ind_max-1)=[];
    end
    ind=find(y~=0);
    y2=y(ind);
    ind2=find(y2>thresh_y);
    y2(ind2)=y2(ind2);
    x2=x(ind)/(tmp1);
    plot(x2,y2,'*-');
%     plot(y2,'*-')
    hold on;
end
xlabel('MDSCN of M component')
ylabel('Number of Coffecients')
legend(img_name{1},img_name{2},img_name{3});

%% N的MDSCN特征的直方图统计
figure('name','N的MDSCN的直方图统计');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    structdis=structdis_N{i};
    ind1=find(structdis>trunc_val*0.05);
    ind2=find(structdis<-trunc_val*0.05);
    ind=[ind1;ind2];
    structdis(ind)=[];
    h=histogram(structdis);
    h_values{i}=h.Values;
    h_edges{i}=h.BinEdges;
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});

% N的折线图
figure('name','N的MDSCN的折线图');
for i=1:cnt_img
    y=h_values{i};
    x=h_edges{i};
    ind=find(y~=0);
    y2=y(ind);
    x2=x(ind)/(tmp1);
    plot(x2,y2,'*-');
    hold on;
end
xlabel('MDSCN of N component')
ylabel('Number of Coffecients')
legend(img_name{1},img_name{2},img_name{3});

%% Lmax的MDSCN特征的直方图统计
figure('name','Lmax的MDSCN的直方图统计');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    structdis=structdis_Lmax{i};
    ind1=find(structdis>trunc_val);
    ind2=find(structdis<-trunc_val);
    ind=[ind1;ind2];
    structdis(ind)=[];
    h=histogram(structdis);
    h_values{i}=h.Values;
    h_edges{i}=h.BinEdges;
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});

% Lmax的折线图
figure('name','Lmax的MDSCN的折线图');
for i=1:cnt_img
    y=h_values{i};
    x=h_edges{i};
    ind=find(y~=0);
    y2=y(ind);
    x2=x(ind)/(tmp1);
    plot(x2,y2,'*-');
    hold on;
end
xlabel('MDSCN of Maximum Local Component')
ylabel('Number of Coffecients')
legend(img_name{1},img_name{2},img_name{3});

%% Lmin的MDSCN特征的直方图统计
figure('name','Lmin的MDSCN的直方图统计');
for i=1:cnt_img
%     subplot(1,cnt_img,i);
    structdis=structdis_M{i};
    ind1=find(structdis>trunc_val*0.01);
    ind2=find(structdis<-trunc_val*0.01);
    ind=[ind1;ind2];
    structdis(ind)=[];
    h=histogram(structdis);
    h_values{i}=h.Values;
    h_edges{i}=h.BinEdges;
    hold on;
end
legend(img_name{1},img_name{2},img_name{3});

% Lmin的折线图
figure('name','Lmin的MDSCN的折线图');
for i=1:cnt_img
    y=h_values{i};
    x=h_edges{i};
    ind=find(y~=0);
    y2=y(ind);
    x2=x(ind)/(tmp1);
    plot(x2,y2,'*-');
    hold on;
end
xlabel('MDSCN of Minimum Local Component')
ylabel('Number of Coffecients')
legend(img_name{1},img_name{2},img_name{3});

%% 显著性特征的MSCN
% figure('name','显著性特征的MSCN');
% for i=1:cnt_img
% mscn_VSMap{i}=MSCN(VSMap{i});
% ind1=find(mscn_VSMap{i}>1);
% ind2=find(mscn_VSMap{i}<-1);
% mscn_VSMap{i}(ind1)=1;
% mscn_VSMap{i}(ind2)=-1;
% subplot(1,cnt_img,i);
% h_mscn_vs{i}=histogram(mscn_VSMap{i});
% legend(img_name{i});
% hold on;
% end


% figure('name','显著性特征MSCN拟合')
% for i=1:cnt_img
%     h_values=h_mscn_vs{i}.Values;
%     h_values2=h_values(1:4:end);
%     plot(h_values2,'*-');
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
%% 原图的直方图
% figure('name','原图的直方图');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(img{i}(:));
%     xlabel(img_name{i});
% end
%% M分量的直方图
% figure('name','M分量的直方图');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(M{i}(:));
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});



%% N分量的直方图
% figure('name','N分量的直方图');
% for i=1:cnt_img
% %     subplot(1,cnt_img,i);
%     histogram(N{i}(:));
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
%% M分量的MSCN
% figure('name','M分量的MSCN');
% for i=1:cnt_img
% mscn_M{i}=MSCN(M{i});
% % subplot(1,cnt_img,i);
% h_mscn_M{i}=histogram(mscn_M{i});
% hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
% 
% figure('name','M分量的MSCN拟合')
% for i=1:cnt_img
%     h_values=h_mscn_M{i}.Values;
%     h_values=medfilt1(h_values,3);
%     plot(h_values,'-');
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
%% N分量的MSCN
% figure('name','N分量的MSCN');
% for i=1:cnt_img
% mscn_N{i}=MSCN(M{i});
% % subplot(1,cnt_img,i);
% h_mscn_N{i}=histogram(mscn_N{i});
% hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});
% 
% figure('name','N分量的MSCN拟合')
% for i=1:cnt_img
%     h_values=h_mscn_N{i}.Values;
%     plot(h_values,'-');
%     hold on;
% end
% legend(img_name{1},img_name{2},img_name{3});

%% 画MSCN的直方图
% figure('name','画MSCN的直方图');
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