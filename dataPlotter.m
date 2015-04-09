datadir = 'bin_debug'; % Location of data files
caseName = 'simData'; % Name of files to be plotted
Iters = [0:200:10000];    % Iterations to be plotted
N = 50;                % # of points to use in each direction for surf plot
GIF = false;           % Save plots to animated GIF file?

data = cell(length(Iters));
for i=1:length(Iters)
    iter = Iters(i);
    filename = sprintf('%s/%s_%09d.vtk',datadir,caseName,iter);
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
    surf(X,Y,U);    zlim([.5,1.1]);
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