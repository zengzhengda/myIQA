% 提取CID2013数据集的特征
function [my_mat]=fetchFeatureAll()
pathd=cell(1,6);
pathd(1)={'..\Datasets\CID2013\IS1\'};
pathd(2)={'..\Datasets\CID2013\IS2\'};
pathd(3)={'..\Datasets\CID2013\IS3\'};
pathd(4)={'..\Datasets\CID2013\IS4\'};
pathd(5)={'..\Datasets\CID2013\IS5\'};
pathd(6)={'..\Datasets\CID2013\IS6\'};
filetype='*.jpg';
% fid=fopen('..\Datasets\CID2013\img_name.txt','w');

fea_mat=[];
for i=1:6
    i
    pad=char(pathd(i));
    dimg_path_list1=dir(pad);
    len_list1=length(dimg_path_list1);
    % dimg_path_list = dir(strcat(pad,filetype));
    for ii=1:len_list1 
        if(dimg_path_list1(ii).name(1)=='c')
            pad2=strcat(pad,dimg_path_list1(ii).name,'\');
            dimg_path_list2=dir(strcat(pad2,filetype));
            len_list2=length(dimg_path_list2);
            for j=1:len_list2
                disname=dimg_path_list2(j).name;       
                %% 记录图像id
%                 if(j~=len_list2)
%                     fprintf(fid,'%s  ',disname);
%                 else
%                     fprintf(fid,'%s\n',disname);
%                 end
                %%
                disname=strcat(pad2,disname);
                disimg=imread(disname);
                fea=fetchFeature(disimg);
                fea_mat=[fea_mat;fea];               

            end
         end
    end
end
% fclose(fid);
path_excel='..\Datasets\CID2013\CID2013 data - version 12112014.xlsx';
[number,txt,raw]=xlsread(path_excel,'CID2013 MOS');
img_mos=number(:,2);
my_mat=[fea_mat,img_mos]; % 最后一列为y
save my_mat.mat my_mat;   
