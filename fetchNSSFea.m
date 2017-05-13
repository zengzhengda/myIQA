function [nss_fea structdis]=fetchNSSFea(disimg)
if(ndims(disimg)==3)
    disimg=rgb2gray(disimg);
end
disimg=double(disimg);
% % 平均滤波
% window=fspecial('gaussian',7,7/6);
% window=window/sum(sum(window));
% mu=filter2(window,disimg,'same');
% mu_sq=mu.*mu;
% sigma=sqrt(abs(filter2(window, disimg.*disimg, 'same') - mu_sq)); % sigma的求法
% structdis=(disimg-mu)./(sigma+eps);
mode =2; % 1 mean,2 median, 3 mean and median 
nss_fea=[];
if(mode==1)
    % 平均滤波
    window=fspecial('gaussian',7,7/6);
    window=window/sum(sum(window));
    mu=filter2(window,disimg,'same');
    mu_sq=mu.*mu;
    sigma=sqrt(abs(filter2(window, disimg.*disimg, 'same') - mu_sq)); % sigma的求法
    structdis=(disimg-mu)./(sigma+eps);
    
    [nss_alpha nss_overallstd]=estimateggdparam(structdis(:));
    nss_skewness=skewness(structdis(:));% 斜度特征
    nss_kurtosis=kurtosis(structdis(:)); % 峰度特征
    nss_entropy=entropy(structdis);% 信息熵特征
    nss_fea=[nss_alpha nss_overallstd nss_skewness nss_kurtosis nss_entropy];
elseif(mode==2)
%采用中值局部归一化
    win_scale=[3,3];
    med_map=medfilt2(disimg,win_scale);
    
    [mad_map]=immovmad(disimg,med_map);
    structdis=(disimg-med_map)./(mad_map+eps);
%     med_sq=med_map.*med_map;
%     tmp1=medfilt2(disimg.*disimg,win_scale);
%     tmp2=tmp1 - med_sq;
%     sigma_rob=sqrt(abs(tmp2));
%     structdis=(disimg-med_map)./(sigma_rob+eps);
%     structdis2=(disimg-med_map);

    [nss_alpha nss_overallstd]=estimateggdparam(structdis(:));
    nss_skewness=skewness(structdis(:));% 斜度特征
    nss_kurtosis=kurtosis(structdis(:)); % 峰度特征
    nss_entropy=entropy(structdis);% 信息熵特征
    nss_fea=[nss_alpha nss_overallstd nss_skewness nss_kurtosis nss_entropy];
else
     % 平均滤波
    window=fspecial('gaussian',7,7/6);
    window=window/sum(sum(window));
    mu=filter2(window,disimg,'same');
    mu_sq=mu.*mu;
    sigma=sqrt(abs(filter2(window, disimg.*disimg, 'same') - mu_sq)); % sigma的求法
    structdis_mean=(disimg-mu)./(sigma+eps);
    
    [nss_alpha_mean nss_overallstd_mean]=estimateggdparam(structdis_mean(:));
    nss_skewness_mean=skewness(structdis_mean(:));% 斜度特征
    nss_kurtosis_mean=kurtosis(structdis_mean(:)); % 峰度特征
    nss_entropy_mean=entropy(structdis_mean);% 信息熵特征
    nss_fea_mean=[nss_alpha_mean nss_overallstd_mean nss_skewness_mean nss_kurtosis_mean nss_entropy_mean];
    
    %采用中值局部归一化
    win_scale=[3,3];
    med_map=medfilt2(disimg,win_scale);
    med_sq=med_map.*med_map;
    sigma_rob=sqrt(abs(medfilt2(disimg.*disimg,win_scale) - med_sq));
    structdis_med=(disimg-med_map)./(sigma_rob+eps);

    [nss_alpha_med nss_overallstd_med]=estimateggdparam(structdis_med(:));
    nss_skewness_med=skewness(structdis_med(:));% 斜度特征
    nss_kurtosis_med=kurtosis(structdis_med(:)); % 峰度特征
    nss_entropy_med=entropy(structdis_med);% 信息熵特征
    nss_fea_med=[nss_alpha_med nss_overallstd_med nss_skewness_med nss_kurtosis_med nss_entropy_med];
    nss_fea=[nss_fea_mean nss_fea_med];
end
end