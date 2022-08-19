function C7_generate_bounds

data = dlmread('bootstrap_simulations.txt');
size(data)
data(data(:,end) == 0,:) = [];
size(data)



dir = '../input/';


cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;


%% read data

bestmodel_lab = dlmread('../output/lab_freq.txt');
bestmodel_neg = dlmread('../output/cluster_size.txt');
bestmodel_lab_hsc = dlmread('../output/real_lab.txt');
bestmodel_neg_hsc = dlmread('../output/real_neg.txt');


data_lab_rel = dlmread([dir 'input_lab_rel_hsc_pv.txt']);
data_lab_number = dlmread([dir 'input_lab_num_tot_nodivide.txt']);
data_cluster_size = dlmread([dir 'input_cluster_size_rel_hsc_over_time_pv.txt']);
data_neg_number = dlmread([dir 'input_neg_num_tot_nodivide.txt']);



%% prepare measured vector


time = data_lab_rel(:,1);

measured_lab_rel = data_lab_rel(:,2:n_clu +1);
errors_lab_rel = data_lab_rel(:,n_clu+2:end);


measured_cluster_size = data_cluster_size(:, 2:n_clu +1);
errors_cluster_size = data_cluster_size(:, n_clu + 2 : end);



measured_lab_num_hsc = data_lab_number(:,1+clu_hsc);
errors_lab_num_hsc = data_lab_number(:,1+clu_hsc+n_clu);
measured_neg_num_hsc = data_neg_number(:,1+clu_hsc);
errors_neg_num_hsc = data_neg_number(:,1+clu_hsc+n_clu);


% t_plot = time(1):270;
t_plot = time(1):0.25:270;


%% read M

M = create_differentiation_matrix(cluster_names,n_clu);


%% generate model for each simulation

ld = size(data,1);

for i = 1:ld
    i
    theta = data(i,1:104);
    
    theta = theta(:);

    
    I_neg = theta(1: n_clu+2);
        
        
        l0 = theta(n_clu+3:(n_clu+1)*2 +2).*I_neg;
        
        d = zeros(n_clu+2,n_clu+2);% differentiation
        
        
        for index = 1: size(M,1)
            
            d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
            
        end
        
        
        
        p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);
        
        p = [p(1:20);sum(d(21,:));p(21:end)];% proliferation
        
        
        
        r = theta(end-1); % logistic parameter
        K = theta(end); % carrying capacity
        
        
        
        k = sum(d,2)-p;
    
    
    % solv lab
    
    sol_lab = ode45(@ODE, [time(1), 270], l0);
    
    
    
    model_lab = deval(sol_lab,t_plot)';
    
    model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu + 1) + model_lab(:,n_clu + 2);
    
    
    
    % solv neg
    
    sol_neg = ode45(@ODE, [time(1), 270], I_neg);
    
    
    model_neg = deval(sol_neg,t_plot)';
    
    model_neg_hsc = model_neg(:,clu_hsc) + model_neg(:,n_clu + 1) + model_neg(:,n_clu + 2);
    
    %
    
    
    
    model_lab(:,clu_hsc) = model_lab_hsc;
    model_neg(:,clu_hsc) = model_neg_hsc;
    
    
    
    %
    
    model_lab(:,n_clu+1:end) = [];
    
    model_neg(:,n_clu+1:end) = [];
    
    
    
    %
    
    model_lab = model_lab ./ repmat(model_lab(:,clu_hsc),1,n_clu);
    model_lab(:,clu_hsc) = [];
    
    model_lab_store(i,:,:) = model_lab;
    
    %
    model_neg = model_neg ./ repmat(model_neg(:,clu_hsc),1,n_clu);
    model_neg(:,clu_hsc) = [];
    
    model_neg_store(i,:,:) = model_neg;
    %
    
    
    
    
    model_lab_hsc_store(i,:) = model_lab_hsc;
    %
    
    model_neg_hsc_store(i,:) = model_neg_hsc;
    
end
size(model_lab_hsc_store)
model_lab_down = quantile(model_lab_store,0.025,1);
model_lab_up = quantile(model_lab_store,0.975,1);


model_neg_down = quantile(model_neg_store,0.025,1);
model_neg_up = quantile(model_neg_store,0.975,1);

model_lab_h_down = quantile(model_lab_hsc_store,0.025,1);
model_lab_h_up = quantile(model_lab_hsc_store,0.975,1);


model_neg_h_down = quantile(model_neg_hsc_store,0.025,1);
model_neg_h_up = quantile(model_neg_hsc_store,0.975,1);


%% plot

figure(9)
clf
for ii = 1:19
    
    subplot(5,4,ii)
    hold on
    
    
    model_up = reshape(model_lab_up(1,:,ii),size(model_lab_up,2),1);
    model_down = reshape(model_lab_down(1,:,ii),size(model_lab_down,2),1);
    
    X = [t_plot';flipud(t_plot')]';
    Y = [model_up;flipud(model_down)]';
    
    
    
    
    
    fill(X,Y,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
    
    
    
    %     plot(t_plot, model_up, 'r', 'Linewidth', 2)
    %     plot(t_plot, model_down, 'b', 'Linewidth', 2)
    %
    plot(bestmodel_lab(:,1),bestmodel_lab(:,ii+2), 'r', 'Linewidth', 1.2)
    
    errorbar(time, measured_lab_rel(:,ii+1), errors_lab_rel(:,ii+1), 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3')
    
    %     line([3,300],[1 1],'color','k')
    
    if ii >16 || ii == 16
        xlabel('time (days)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('Rel. label frequency')
    end
    
    title(cluster_names(ii + 1))
    
    
    set(gca,'xscale','log')
    set(gca,'XTickLabel',[1, 10, 100])

end



set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/lab_freq_bootstrap.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width', 22,'height', 25,'Renderer','painters','Lockaxes',0);%

%%



figure(98)
clf

hold on

for ii = 1:19
    
    subplot(5,4,ii)
    hold on
    
    
    model_up = reshape(model_neg_up(1,:,ii),size(model_neg_up,2),1);
    model_down = reshape(model_neg_down(1,:,ii),size(model_neg_down,2),1);
    
    X = [t_plot';flipud(t_plot')]';
    Y = [model_up;flipud(model_down)]';
    
    
    
    fill(X,Y,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
    
    
    plot(bestmodel_neg(:,1),bestmodel_neg(:,ii+2), 'r', 'Linewidth', 1.2)
    
    
    
    errorbar(time, measured_cluster_size(:,ii+1), errors_cluster_size(:,ii+1), 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)
    
    
    if ii >16 || ii == 16
        xlabel('time (days)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('Rel. cluster size')
    end
    
    title(cluster_names(ii+1))
    
    
    set(gca,'xscale','log')
    set(gca,'XTickLabel',[1,10, 100])

end

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/cluster_size_bootstrap.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width',22,'height',25,'Renderer','painters','Lockaxes',0);%


%%

figure(1001)
clf

subplot(1,2,1)

hold on


  model_up = model_lab_h_up';
    model_down = model_lab_h_down';
    
    X = [t_plot';flipud(t_plot')]';
    Y = [model_up;flipud(model_down)]';
    
 
    
    fill(X,Y,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
    



    plot(bestmodel_lab_hsc(:,1),bestmodel_lab_hsc(:,2), 'r', 'Linewidth', 1)

errorbar(time,measured_lab_num_hsc,errors_lab_num_hsc, 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)




xlabel('time (days)')
ylabel ('No of labelled cells')

% legend('data','model')

%

subplot(1,2,2)

hold on


  model_up = model_neg_h_up';
    model_down = model_neg_h_down';
    
    X = [t_plot';flipud(t_plot')]';
    Y = [model_up;flipud(model_down)]';
    
 
    
    fill(X,Y,'b','facecolor',[1 1 1]*0.8,'edgecolor',[1 1 1]*0.8)
    

    plot(bestmodel_neg_hsc(:,1),bestmodel_neg_hsc(:,2), 'r', 'Linewidth', 1)

errorbar(time,measured_neg_num_hsc,errors_neg_num_hsc, 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)




xlabel('time (days)')
ylabel('Cluster size')


% legend('data','model')
set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/real_lab_neg_bootstrap.eps','FontMode', 'fixed','Fontsize',6 ,'color', 'cmyk','width',13,'height',4,'Renderer','painters','Lockaxes',0);%


%%
    function dxdt = ODE(~,x)
        
        
         
            dxdt = d' * x  - k .* x;
            
            xx = x(1)+x(21)+x(22);
            
            dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));        
        
        
        
    end

end