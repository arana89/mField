%% example of using mField class

obj = mField();
%% show available methods
methods(obj)
%% show available properties
properties(obj)
%% use loadOVF to get magnetization field of Hopfion and assign to mx,my,mz
obj.loadOVF('./testData/hopfion_n25_r004_H1_z0.ovf');
%%
obj.showPlot = 1; %boolean to plot Hopfion
obj.hopfNumber;
%% print computed Hopf number
obj.tHopfion