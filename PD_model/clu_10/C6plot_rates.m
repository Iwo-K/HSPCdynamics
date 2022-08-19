load('./PD_matrixes/data_input_pd_clu_10.mat')
load('./PD_matrixes/results_clu_10.mat')


solbest = 2;


par = computeParameters(Est{solbest}.parameters.MS.par(:,1));
%% compute drift
n_grid1 = 300;
grid = linspace(0,1,n_grid1);
par.grid = grid(2:end-1);
parameters = Est{solbest}.parameters.MS.par(:,1);

grid_x = linspace(grid(1)+(grid(2)-grid(1))/2,grid(end-1)+...
    (grid(end)-grid(end-1))/2,n_grid1-1); 
par.gridx = grid_x;

v = parameters(10:18);

x = (0:0.125:1)';

y = v(:);

xx = grid(2:end-1);

yy = sym(size(xx));

h = x(2)-x(1);
b = 1/h*(y(2:end)-y(1:end-1));
v = 4*h;
u = 6*(b(2:end)-b(1:end-1));


z = [0;inv(diag(v*ones(length(x)-2,1))+diag(h*ones(length(x)-3,1),-1)+diag(h*ones(length(x)-3,1),1))*u;0];



for i=1:length(xx)
    ind = min(find(x<=xx(i),1,'last'),length(x)-1);
    yy(i) = z(ind+1)/(6*h)*(xx(i)-x(ind))^3+z(ind)/(6*h)*(x(ind+1)-xx(i))^3+...
        (y(ind+1)/h-z(ind+1)/6*h)*(xx(i)-x(ind))+(y(ind)/h-h/6*z(ind))*(x(ind+1)-xx(i));
   
end

drift_no_corr = double(exp(yy));
%% plot parameters splines


figure(10020)


clf

subplot(1,3,1)
hold on
plot(par.grid,par.D)
xlabel('cell state')
ylabel('diffusion parameter (d^{-1})')

subplot(1,3,2)
hold on

plot(par.grid,drift_no_corr)
xlabel('cell state')
ylabel('drift parameter (d^{-1})')


subplot(1,3,3)
hold on

plot(par.gridx,par.a)
xlabel('cell state')
ylabel('birth-death parameter (d^{-1})')


set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/rates_clu_10.eps','FontMode', 'fixed','Fontsize',13,'color', 'cmyk','width',30,'height',10,'Renderer','painters','Lockaxes',0);%


mat_parameters = [par.grid;par.D;drift_no_corr]';
mat_parameters2 = [par.gridx; par.a]';

dlmwrite('./tables/parameters_clu_10_dpt_diff_drift.txt',mat_parameters)
dlmwrite('./tables/parameters_clu_10_dpt_growth.txt',mat_parameters2)
