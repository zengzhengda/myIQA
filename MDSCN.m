%Í¼Ïñ±ê×¼»¯
function [mdscn_map]=MDSCN(disimg)
if(ndims(disimg)==3)
    disimg=rgb2gray(disimg);
end
disimg=double(disimg);
win_scale=[3,3];
med_map=medfilt2(disimg,win_scale);
med_sq=med_map.*med_map;
sigma_rob=sqrt(abs(medfilt2(disimg.*disimg,win_scale) - med_sq));
mdscn_map=(disimg-med_map)./(sigma_rob+eps);