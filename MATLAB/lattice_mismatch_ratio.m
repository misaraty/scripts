clc;clear;
[filepath]=fileparts([mfilename('fullpath'),'.m']); 
addpath(genpath(filepath))
cd(filepath)

a=1:10;b=1:10;
[a,b]=meshgrid(a,b);
y1=(5.41693.*a-3.9102.*b)./(5.41693.*a);
y2=(5.41693.*a-3.9102.*b)./(3.9102.*b);
c=find(abs(y1)<=0.1&abs(y2)<=0.1);
[a(c) b(c)]

