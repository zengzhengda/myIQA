function [nss_alpha nss_overallstd nss_skewness nss_kurtosis nss_entropy structdis]=fetchNSSFea(disimg)
if(ndims(disimg)==3)
    disimg=rgb2gray(disimg);
end
disimg=double(disimg);
% % ƽ���˲�
% window=fspecial('gaussian',7,7/6);
% window=window/sum(sum(window));
% mu=filter2(window,disimg,'same');
% mu_sq=mu.*mu;
% sigma=sqrt(abs(filter2(window, disimg.*disimg, 'same') - mu_sq)); % sigma����
% structdis=(disimg-mu)./(sigma+eps);

%������ֵ�ֲ���һ��
win_scale=[4,4];
med_map=medfilt2(disimg,win_scale);
med_sq=med_map.*med_map;
sigma_rob=sqrt(abs(medfilt2(disimg.*disimg,win_scale) - med_sq));
structdis=(disimg-med_map)./(sigma_rob+eps);

[nss_alpha nss_overallstd]=estimateggdparam(structdis(:));

nss_skewness=skewness(structdis(:));% б������
nss_kurtosis=kurtosis(structdis(:)); % �������
nss_entropy=entropy(structdis);% ��Ϣ������
% 
% % б�Ⱥͷ������
% shifts                   = [ 0 1;1 0 ; 1 1; -1 1];
%  
% for itr_shift =1:4
% shifted_structdis        = circshift(structdis,shifts(itr_shift,:));
% pair                     = structdis(:).*shifted_structdis(:);
% nss_skewness1=skewness(pair);% б������
% nss_kurtosis1=kurtosis(pair);% �������
% nss_skewness=[nss_skewness nss_skewness1];
% nss_kurtosis=[nss_kurtosis nss_kurtosis1];
% 
% % [alpha leftstd rightstd] = estimateaggdparam(pair);
% % const                    =(sqrt(gamma(1/alpha))/sqrt(gamma(3/alpha)));
% % meanparam                =(rightstd-leftstd)*(gamma(2/alpha)/gamma(1/alpha))*const;
% % feat                     =[feat alpha meanparam leftstd^2 rightstd^2];
% end
%     
end