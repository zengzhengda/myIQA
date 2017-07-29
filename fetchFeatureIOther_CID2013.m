function []=fetchFeatureOther_CID2013(metric)
%% 提取CID2013数据集的其他方法的特征
	pathd=cell(1,6);
	pathd(1)={'..\Datasets\CID2013\IS1\'};
	pathd(2)={'..\Datasets\CID2013\IS2\'};
	pathd(3)={'..\Datasets\CID2013\IS3\'};
	pathd(4)={'..\Datasets\CID2013\IS4\'};
	pathd(5)={'..\Datasets\CID2013\IS5\'};
	pathd(6)={'..\Datasets\CID2013\IS6\'};
	filetype='*.jpg';
	% fid=fopen('..\Datasets\CID2013\img_name.txt','w');

	len_imgs=474;
	ori_img_dataset=cell(1,len_imgs);
	fea_mat=[];
	cnt_img=1; % 图像计数
	% 主观数据集
	path_excel='..\Datasets\CID2013\CID2013 data - version 12112014.xlsx';
	[number,txt,raw]=xlsread(path_excel,'CID2013 MOS');
	img_mos=number(:,3);

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
	                disname_all=strcat(pad2,disname);
	                disimg=imread(disname_all);   
	                %% 重新保存整理图像
	%                 img_gray=rgb2gray(disimg);
	%                 path_gray=strcat('..\Datasets\CID2013_gray\',disname(1:end-4),'_',num2str(img_mos(cnt_img)),'.jpg');
	%                 path_color=strcat('..\Datasets\CID2013_color\',disname(1:end-4),'_',num2str(img_mos(cnt_img)),'.jpg');
	%                 imwrite(img_gray,path_gray);
	%                 imwrite(disimg,path_color);
	                %%
	                ori_img_dataset{cnt_img}=disimg;
	                cnt_img=cnt_img+1;
	                switch metric
	                	case 'IDEAL'
	                		fea=fetchFeatureIDEAL(disimg); % 分析salient
	                		fea_mat=[fea_mat;fea];
	                	otherwise
	                		print('error');
	                		break;
	                end
	                               

	            end
	         end
	    end
	end
	% fclose(fid);
	img_ind=(1:len_imgs)';
	my_mat_ideal=[img_ind,fea_mat,img_mos]; % 最后一列为y
	save ideal_mat_cid2013.mat my_mat_ideal;   
	% save ori_img_dataset.mat ori_img_dataset;
end