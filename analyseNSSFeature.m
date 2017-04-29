% 2017-04-25 修改
function [feature]=analyseNSSFeature(disimg,cnt_img_level)
[rows,cols,dim]=size(disimg);
if(nargin==1)
    cnt_img_level=1;
end
feature=[];
for level=1:cnt_img_level

   %% 图像灰度化
    if(size(disimg,3)~=1)
        disimg_gray=double(rgb2gray(disimg));
    end
    %% 颜色空间转化：RGB->LMN 基本特征
    L = 0.06 * double(disimg(:,:,1)) + 0.63 * double(disimg(:,:,2)) + 0.27 * double(disimg(:,:,3));
    M = 0.30 * double(disimg(:,:,1)) + 0.04 * double(disimg(:,:,2)) - 0.35 * double(disimg(:,:,3));
    N = 0.34 * double(disimg(:,:,1)) - 0.60 * double(disimg(:,:,2)) + 0.17 * double(disimg(:,:,3));
%     histogram(L);
    %% MN分量的统计特征
    [M_nss_fea]=fetchNSSFea(M);
    [N_nss_fea]=fetchNSSFea(N);
    
    %% NSS特征
    [L_nss_fea]=fetchNSSFea(disimg_gray);
    %% 特征向量
    % MN颜色特征10项 NSS特征5
    feature=[feature  L_nss_fea M_nss_fea N_nss_fea];
    
    disimg=imresize(disimg,0.5);
end