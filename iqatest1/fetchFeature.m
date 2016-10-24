function [feature]=fetchFeature(disimg)
[rows,cols,dim]=size(disimg);
L = 0.06 * double(disimg(:,:,1)) + 0.63 * double(disimg(:,:,2)) + 0.27 * double(disimg(:,:,3));
M = 0.30 * double(disimg(:,:,1)) + 0.04 * double(disimg(:,:,2)) - 0.35 * double(disimg(:,:,3));
N = 0.34 * double(disimg(:,:,1)) - 0.60 * double(disimg(:,:,2)) + 0.17 * double(disimg(:,:,3));

%% LMN 的期望和方差
mu_L=mean(L(:));
temp=(L-mu_L).^2;
sigma_L=mean(temp(:));
mu_M=mean(M(:));
temp=(M-mu_M).^2;
sigma_M=mean(temp(:));
mu_N=mean(N(:));
temp=(N-mu_N).^2;
sigma_N=mean(temp(:));
feature1=[mu_L,sigma_L,mu_M,sigma_M,mu_N,sigma_N];

%% FISH特征
if(size(disimg,3)~=1)
    disimg=rgb2gray(disimg);
end
fish_fea=fetch_fish_fea(disimg);

%% 特征向量
feature=[fish_fea feature1];