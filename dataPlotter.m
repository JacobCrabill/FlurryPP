clear all; close all;

datadir = 'bin'; % Location of data files
caseName = 'simData'; % Name of files to be plotted
Iters = [0:500:50000];    % Iterations to be plotted
N = 50;                % # of points to use in each direction for surf plot
GIF = false;           % Save plots to animated GIF file?

data = cell(length(Iters));
for i=1:length(Iters)
    iter = Iters(i);
    filename = sprintf('%s/%s.csv.%09d',datadir,caseName,iter);
    data{i} = csvread(filename,1,0);
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
    u = data{i}(:,6); % x y z rho u v p
    F = TriScatteredInterp(x,y,u);
    
    U = log(abs(F(X,Y)-10));
    
    %surf(X,Y,U);    zlim([.5,1.1]);
    %contourf(X,Y,U,100,'edgecolor','none'); colorbar;
    pcolor(X,Y,U); shading flat; colorbar;
    pause(.1);
    %pause;
    
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