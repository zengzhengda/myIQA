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
    [M_nss_fea]=fetchNSSFea(M);
    [N_nss_fea]=fetchNSSFea(N);
    
    %% NSS����
    [L_nss_fea]=fetchNSSFea(disimg_gray);
    %% ��������
    % MN��ɫ����10�� NSS����5
    feature=[feature  L_nss_fea M_nss_fea N_nss_fea];
    
    disimg=imresize(disimg,0.5);
end