% 2017-04-25 �޸�
function [feature]=analyseNSSFeature(disimg,cnt_img_level)
[rows,cols,dim]=size(disimg);
if(nargin==1)
    cnt_img_level=1;
end
feature=[];
for level=1:cnt_img_level

   %% ͼ��ҶȻ�
    if(size(disimg,3)~=1)
        disimg_gray=double(rgb2gray(disimg));
    end
    %% ��ɫ�ռ�ת����RGB->LMN ��������
    L = 0.06 * double(disimg(:,:,1)) + 0.63 * double(disimg(:,:,2)) + 0.27 * double(disimg(:,:,3));
    M = 0.30 * double(disimg(:,:,1)) + 0.04 * double(disimg(:,:,2)) - 0.35 * double(disimg(:,:,3));
    N = 0.34 * double(disimg(:,:,1)) - 0.60 * double(disimg(:,:,2)) + 0.17 * double(disimg(:,:,3));
%     histogram(L);
    %% MN������ͳ������
    [M_alpha M_overallstd M_skewness M_kurtosis M_entropy]=fetchNSSFea(M);
    [N_alpha N_overallstd N_skewness N_kurtosis N_entropy]=fetchNSSFea(N);
    fea_MN=[M_alpha M_overallstd M_skewness M_kurtosis M_entropy N_alpha N_overallstd N_skewness N_kurtosis N_entropy];
    
    %% NSS����
    [nss_alpha ,nss_overallstd, nss_skewness, nss_kurtosis ,nss_entropy]=fetchNSSFea(disimg_gray);
    feature_nss=[nss_alpha nss_overallstd nss_skewness nss_kurtosis nss_entropy]; % 5��
    %% ��������
    % MN��ɫ����10�� NSS����5
    feature=[feature fea_MN feature_nss];
    
    disimg=imresize(disimg,0.5);
end