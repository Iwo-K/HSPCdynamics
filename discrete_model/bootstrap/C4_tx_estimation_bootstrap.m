function C4_tx_estimation_bootstrap
%% bootstrap results

data_boot = dlmread('bootstrap_simulations.txt');

data_boot(data_boot(:,end) == 0,:) = [];

time = 3:0.1:8;

%%
input_dir = '../../procdata/07script/';
cluster_names = [0:12,14,16,20,24,25,26,28];
n_clu = length(cluster_names);
clu_hsc = 1;

%Transplantation data - counts per cluster
data = readtable([input_dir 'tx_cellnos.csv']);
data_m = [data.x1Day, data.x3Day, data.x5Day, data.x7Day]';
data_norm = data_m ./ repmat(data_m(:,clu_hsc),1,n_clu);

err_data_norm = dlmread('./input_tx_boot_err_prop.txt');


% Loading cluster sizes
best = dlmread('../output/best.txt');
best(end) = [];
clu_size = best(1:22);
hsc_totsize = sum(clu_size([clu_hsc n_clu+1 n_clu+2]));
cl0_relsize = clu_size(clu_hsc) / hsc_totsize;
cl30_relsize = clu_size(n_clu + 1) / hsc_totsize;
cl40_relsize = clu_size(n_clu + 2) / hsc_totsize;
%% read M

M = create_differentiation_matrix(cluster_names,n_clu);

%% define parameters
    theta = best(:);
    
    d = zeros(n_clu+2,n_clu+2);
    
    
    for index = 1: size(M,1)
        
        d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
        
    end
    
    
    p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);
    
    p = [p(1:20);sum(d(21,:));p(21:end)];

    
    r = theta(end-1);
    K = theta(end);
    
    
    
    k = sum(d,2)-p;

%% Day3 as starting point

c3 = data.x3Day;
c3 = c3/c3(1);
% Splitting cluster 0 into: 0, 30, 40 with in the proportions from the Hoxb5
% model
c3(clu_hsc) = cl0_relsize;
c3 = [c3; cl30_relsize; cl40_relsize];

% time = [1 3 5 7];
sol = ode45(@ODE, [time(1), 8], c3);
model3_best = deval(sol,time)';

% %Re-normalise to cluster 0 (sum of cluster 0,30,40)

hsc = model3_best(:,clu_hsc) + model3_best(:,n_clu+1) + model3_best(:,n_clu+2);
model3_best(:,clu_hsc) = hsc;
model3_best(:,n_clu+1:end) = [];
model3_best = model3_best ./ repmat(model3_best(:,clu_hsc),1,n_clu);

%%
ld = size(data_boot,1);

for i = 1:ld
    
    i
    
    theta = data_boot(i,1:104);
    
    theta = theta(:);
    
    d = zeros(n_clu+2,n_clu+2);
    
    
    for index = 1: size(M,1)
        
        d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
        
    end
    
    
    p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);
    
    p = [p(1:20);sum(d(21,:));p(21:end)];
    
    
    r = theta(end-1);
    K = theta(end);
    
    
    k = sum(d,2)-p;
    
    
    sol = ode45(@ODE, [time(1), 8], c3); %check?
    model3 = deval(sol,time)';
    
    %Re-normalise to cluster 0 (sum of cluster 0,30,40)
    hsc = model3(:,clu_hsc) + model3(:,n_clu+1) + model3(:,n_clu+2);
    model3(:,clu_hsc) = hsc;
    model3(:,n_clu+1:end) = [];
    model3_store(i,:,:) = model3 ./ repmat(model3(:,clu_hsc),1,n_clu);
    
end

%%

model_down = quantile(model3_store,0.025,1);
model_up = quantile(model3_store,0.975,1);


%% Day3 model figure
figure(2)
clf

hold on

for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    
    
    model_u = reshape(model_up(1,:,ii),size(model_up,2),1);
    model_d = reshape(model_down(1,:,ii),size(model_down,2),1);
    
    X = [time';flipud(time')]';
    Y = [model_u;flipud(model_d)]';
    fill(X,Y,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
    
    
    plot(time, model3_best(:, ii), 'r', 'Linewidth', 1)
    errorbar([3 5 7], data_norm([2,3,4], ii), err_data_norm([2,3,4], ii),  'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)


    if ii >16 || ii == 16
        xlabel('Time (days)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('Rel. cluster size')
    end
    
    title(cluster_names(ii))
    
    
    set(gca,'xscale','log')
end



set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/transplant_bootstrap.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width',22,'height',15,'Renderer','painters','Lockaxes',0);%

%% ODE

    function dxdt = ODE(~,x)
        
        dxdt = d' * x  - k .* x;
        
        xx = x(1)+x(21)+x(22);
        
        dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
        
    end


end
