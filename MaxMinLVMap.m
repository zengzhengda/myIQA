function [map_max map_min] = MaxMinLVMap(img)

[xs, ys] = size(img);
x=img;

x1=zeros(xs,ys);
x2=zeros(xs,ys);
x3=zeros(xs,ys);
x4=zeros(xs,ys);
x5=zeros(xs,ys);
x6=zeros(xs,ys);
x7=zeros(xs,ys);
x8=zeros(xs,ys);
x9=zeros(xs,ys);

x1(1:xs-2,1:ys-2) = x(2:xs-1,2:ys-1);
x2(1:xs-2,2:ys-1) = x(2:xs-1,2:ys-1);
x3(1:xs-2,3:ys)   = x(2:xs-1,2:ys-1);
x4(2:xs-1,1:ys-2) = x(2:xs-1,2:ys-1);
x5(2:xs-1,2:ys-1) = x(2:xs-1,2:ys-1);
x6(2:xs-1,3:ys)   = x(2:xs-1,2:ys-1);
x7(3:xs,1:ys-2)   = x(2:xs-1,2:ys-1);
x8(3:xs,2:ys-1)   = x(2:xs-1,2:ys-1);
x9(3:xs,3:ys)     = x(2:xs-1,2:ys-1);

x1=x1(2:xs-1,2:ys-1);
x2=x2(2:xs-1,2:ys-1);
x3=x3(2:xs-1,2:ys-1);
x4=x4(2:xs-1,2:ys-1);
x5=x5(2:xs-1,2:ys-1);
x6=x6(2:xs-1,2:ys-1);
x7=x7(2:xs-1,2:ys-1);
x8=x8(2:xs-1,2:ys-1);
x9=x9(2:xs-1,2:ys-1);

d1=(x1-x5);
d2=(x2-x5);
d3=(x3-x5);
d4=(x4-x5);
d5=(x6-x5);
d6=(x7-x5);
d7=(x8-x5);
d8=(x9-x5);

% 最大局部图
dd=max(d1,d2);
dd=max(dd,d3);
dd=max(dd,d4);
dd=max(dd,d5);
dd=max(dd,d6);
dd=max(dd,d7);
dd=max(dd,d8);
% 最小局部图
dd2=max(dd2,d2);
dd2=max(dd2,d3);
dd2=max(dd2,d4);
dd2=max(dd2,d5);
dd2=max(dd2,d6);
dd2=max(dd2,d7);
dd2=max(dd2,d8);
map_max = dd;
map_min = dd2;
end % function





