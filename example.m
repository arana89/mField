%% example of using mField class

obj = mField();
%% show available methods
methods(obj)
%% show available properties
properties(obj)
%% use loadOVF to get magnetization field of Hopfion 
obj.loadOVF('./testData/hopfion_n25_r004_H1_z0.ovf');
%%
obj.showPlot = 1; %boolean to plot Hopfion
obj.hopfNumber();
%% plot B field streamlines
nSeed = 50; %optional input parameter
obj.streamB(nSeed)
%% print computed Hopf number
obj.tHopfion
%% assign physical quantities and compute mach number
obj.Ms = 790e3;
obj.Aex = 10e-12;
obj.cellSize = 5e-9 * [1 1 1];

obj.machNumber(); %assigned to obj.tMach

