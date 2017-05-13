function [mad_map]=immovmad(img,med_map,x)
% Í¼ÏñµÄmedian absolute deviation
% med_map is the median map, x is the half long of side window 
if(nargin <3)
    x=1;
end
if(size(img,3)~=1)
    img=rgb2gray(img);
end
[m,n]=size(img);
img_pad=zeros(m+x*2,n+x*2);
img_pad(x+1:x+m,x+1:n+x)=img;
mad_map=zeros(m,n);
for i=x+1:x+m
    for j=x+1:x+n
        patch=img_pad(i-x:i+x,j-x:j+x);
        tmp=abs(patch-med_map(i-x,j-x));
        mad_map(i-x,j-x)=mean2(tmp(:));
    end
end
end