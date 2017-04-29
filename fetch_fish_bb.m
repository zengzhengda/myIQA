function [fish_bb]=fetch_fish_bb(img)
if(ndims(img)==3)
    img=rgb2gray(img);
end
%% ����ֲ�fish
[M,N]=size(img);
size_blk=64; % ��ĳߴ�
gap=32; % ��֮��ļ����С
fish_map=[];
for row=1:gap:(floor(M/gap)-2)*gap
    for col=1:gap:(floor(N/gap)-2)*gap
        img_blk=img(row:row+size_blk,col:col+size_blk);
        fish_loc=fetch_fish(img_blk);
        fish_map=[fish_map,fish_loc];
    end
end
fish_bb=sqrt(mean(fish_map.^2));
end
