close all; clear all;

Iters = [0:10:500];
datadir = 'bin_debug';
N = 50; % # of points to use in each direction for surf plot

data = cell(length(Iters));
for i=1:length(Iters)
    iter = Iters(i);
    filename = sprintf('%s/simData_%09d.vtk',datadir,iter);
    data{i} = csvread(filename);
end

%% Assuming constant range in x,y (if not, move this inside the 'for' loop)
x = data{i}(:,1);
y = data{i}(:,2);
xx = min(x):(max(x)-min(x))/N:max(x);
yy = min(y):(max(y)-min(y))/N:max(y);
[X,Y] = meshgrid(xx,yy);

for i=1:length(Iters)
    x = data{i}(:,1);
    y = data{i}(:,2);
    u = data{i}(:,3);
    F = TriScatteredInterp(x,y,u);
    
    U = F(X,Y);
    surf(X,Y,U);
    pause(.1);
end