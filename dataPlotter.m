% clear all; close all;
% 
% datadir = 'bin'; % Location of data files
% caseName = 'testRun'; % Name of files to be plotted
% Iters = [0:100:10000];    % Iterations to be plotted
% N = 300;                % # of points to use in each direction for surf plot
% GIF = false;           % Save plots to animated GIF file?
% 
% data = cell(length(Iters));
% for i=1:length(Iters)
%     iter = Iters(i);
%     filename = sprintf('%s/%s.csv.%09d',datadir,caseName,iter);
%     data{i} = csvread(filename,1,0);
% end

%% Assuming constant range in x,y (if not, move this inside the 'for' loop)
x = data{1}(:,1);
y = data{1}(:,2);
%xx = min(x):(max(x)-min(x))/N:max(x);
%yy = min(y):(max(y)-min(y))/N:max(y);
xx = -5:(20/N):15;
yy = -10:(20/N):10;
[X,Y] = meshgrid(xx,yy);

for i=1:length(Iters)
    % x y z rho u v p
    x = data{i}(:,1);
    y = data{i}(:,2);
    rho = data{i}(:,4);
    if size(data{i},2) > 4
        u = data{i}(:,5);
        v = data{i}(:,6);
        p = data{i}(:,7);
        M = sqrt(u.^2+v.^2)./sqrt(1.4*p./rho);
        
        % Choose which quantity to plot
        q = rho;
    else
        q = rho;
    end  
    
    % Interpolate Data to Grid
    F = TriScatteredInterp(x,y,q);
    U = F(X,Y);

    pcolor(X,Y,U); shading flat; colorbar;
    pause(.1);
    
	% Create animated GIF
    if GIF
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end