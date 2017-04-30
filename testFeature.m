clear ;
clc;
close all;
warning('off'); % ≤ªœ‘ æwarning

img={};
path{1}='IS_I_C03_D04_7.3705';
path{2}='IS_I_C03_D01_31.0333';
path{3}='IS_I_C03_D02_94.3333';

img_name{1}='IS\_I\_C03\_D04\_7.3705';
img_name{2}='IS\_I\_C03\_D01\_31.0333';
img_name{3}='IS\_I\_C03\_D02\_94.3333';

for i=1:length(path)
    img{i}=imread([path{i},'.jpg']);
end
cnt_img=length(img);
cnt_img_level=1;
for i=1:cnt_img
[fea]=fetchFeature3(img{i},cnt_img_level);
end