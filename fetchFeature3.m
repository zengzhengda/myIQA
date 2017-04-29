% 2017-04-25 修改
function [feature]=fetchFeature3(disimg,cnt_img_level)
[rows,cols,dim]=size(disimg);
if(nargin==1)
    cnt_img_level=1;
end
feature=[];
for level=1:cnt_img_level
     %% 图像预处理
    %  disimg=imgPreDeal(disimg);

   %% 图像灰度化
    if(size(disimg,3)~=1)
        disimg_gray=double(rgb2gray(disimg));
    end
    %% 颜色空间转化：RGB->LMN 基本特征
    L = 0.06 * double(disimg(:,:,1)) + 0.63 * double(disimg(:,:,2)) + 0.27 * double(disimg(:,:,3));
    M = 0.30 * double(disimg(:,:,1)) + 0.04 * double(disimg(:,:,2)) - 0.35 * double(disimg(:,:,3));
    N = 0.34 * double(disimg(:,:,1)) - 0.60 * double(disimg(:,:,2)) + 0.17 * double(disimg(:,:,3));
%     histogram(L);
    
    window = fspecial('gaussian', 7, 1.5);
    window = window/sum(sum(window));  
    mu_L_map=filter2(window,L,'same');
    mu_M_map=filter2(window,M,'same');
    mu_N_map=filter2(window,N,'same');
    % 期望和方差
    mu_L=mean2(mu_L_map);
    mu_M=mean2(mu_M_map);
    mu_N=mean2(mu_N_map);
    std_L=std2(mu_L_map);
    std_M=std2(mu_M_map);
    std_N=std2(mu_N_map);
    feature_base=[mu_L,mu_M,mu_N,std_L,std_M,std_N];% 6项
    % 鲁棒统计量 
    med_L=median(mu_L_map(:));
    med_M=median(mu_M_map(:));
    med_N=median(mu_N_map(:));
    mad_L=mean2(abs(mu_L_map-med_L));
    mad_M=mean2(abs(mu_M_map-med_M));
    mad_N=mean2(abs(mu_N_map-med_N));
    med_lmn=[med_L med_M med_N];
    mad_lmn=[mad_L mad_M mad_N];
    feature_base_robust=[med_lmn mad_lmn]; %6项
     %% 梯度特征
    [ave_grad std_grad med_grad mad_grad grad_map1]=fetchGradient(disimg_gray);
    grad_fea=[ave_grad std_grad]; % 6项
    grad_fea_robust=[med_grad mad_grad]; % 6项
    %% 显著性特征
    sigmaF = 1.34;%fixed donot change
    omega0 = 0.0210;%fixed
    sigmaD = 145;%fixed
    sigmaC = 0.001;%fixed
    [mu_saliency,sigma_saliency,beta_saliency,VSMap]=fetchSaliencyFeature(disimg,sigmaF,omega0,sigmaD,sigmaC);
    feature_saliency=[mu_saliency,beta_saliency];   
    %% MN分量的统计特征
    [M_alpha M_overallstd M_skewness M_kurtosis M_entropy]=fetchNSSFea(M);
    [N_alpha N_overallstd N_skewness N_kurtosis N_entropy]=fetchNSSFea(N);
    fea_MN=[M_alpha M_overallstd M_skewness M_kurtosis M_entropy N_alpha N_overallstd N_skewness N_kurtosis N_entropy];
    
    %% FISH特征
    fish_fea=fetch_fish(disimg_gray);% 1项
    %% MLV特征（描述sharpness）
    mlv_fea=MLVSharpnessMeasure(disimg_gray);
    %% sharpness 特征
    sharpness_fea=[mlv_fea fish_fea]; % 3项
%     sharpness_fea=[mlv_fea grad_fea]; % 8项
    %% NSS特征
    [nss_alpha ,nss_overallstd, nss_skewness, nss_kurtosis ,nss_entropy]=fetchNSSFea(disimg_gray);
    feature_nss=[nss_alpha nss_overallstd nss_skewness nss_kurtosis nss_entropy]; % 5项
%     %% 斜度特征
%     skewness_fea=skewness(disimg_gray(:));
%     %% 峰度特征
%     kurtosis_fea=kurtosis(disimg_gray(:));
%     %% 信息熵特征
%     entropy_fea=entropy(disimg_gray);
    %% 特征向量
    % 基本特征6 鲁棒基本特征6 梯度特征6 鲁棒梯度特征6 sharpness特征3项 MN颜色特征10项 NSS特征5  显著性特征2
    % 共44项
    feature=[feature feature_base grad_fea feature_base_robust grad_fea_robust fea_MN sharpness_fea feature_nss  feature_saliency];
    
    disimg=imresize(disimg,0.5);
end