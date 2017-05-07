% 2017-04-25 �޸�
function [feature_saliency]=analyseSalientFeature(disimg,cnt_img_level)
[rows,cols,dim]=size(disimg);
if(nargin==1)
    cnt_img_level=1;
end
feature_saliency=[];
for level=1:cnt_img_level

   %% ͼ��ҶȻ�
    if(size(disimg,3)~=1)
        disimg_gray=double(rgb2gray(disimg));
    end
    %% ����������
    sigmaF = 1.34;%fixed donot change
    omega0 = 0.0210;%fixed
    sigmaD = 145;%fixed
    sigmaC = 0.001;%fixed
    [mu_saliency,sigma_saliency,beta_saliency,VSMap]=fetchSaliencyFeature(disimg,sigmaF,omega0,sigmaD,sigmaC);
    feature_saliency=[feature_saliency beta_saliency mu_saliency,sigma_saliency];   
    
    disimg=imresize(disimg,0.5);
end