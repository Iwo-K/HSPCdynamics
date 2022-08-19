function C1_main
%% general info

dir = './input/';

name_best = './output/best.txt';

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;

%% read data

data_lab_rel = dlmread([dir 'input_lab_rel_hsc_pv.txt']);
data_lab_number = dlmread([dir 'input_lab_num_tot_nodivide.txt']);
data_cluster_size = dlmread([dir 'input_cluster_size_rel_hsc_over_time_pv.txt']);
data_neg_number = dlmread([dir 'input_neg_num_tot_nodivide.txt']);

%% prepare measured vector

time = data_lab_rel(:,1);

measured_lab_rel_temp = data_lab_rel(:,2:n_clu +1);
measured_lab_rel = measured_lab_rel_temp;
measured_lab_rel(:,clu_hsc) = [];

errors_lab_rel_temp = data_lab_rel(:,n_clu+2:end);
errors_lab_rel = errors_lab_rel_temp;
errors_lab_rel(:,clu_hsc) = [];


measured_cluster_size_temp = data_cluster_size(:, 2:n_clu +1);
measured_cluster_size = measured_cluster_size_temp;
measured_cluster_size(:,clu_hsc) = [];

errors_cluster_size_temp = data_cluster_size(:, n_clu + 2 : end);
errors_cluster_size = errors_cluster_size_temp;
errors_cluster_size(:,clu_hsc) = [];


measured_lab_num_hsc = data_lab_number(:,1+clu_hsc);
errors_lab_num_hsc = data_lab_number(:,1+clu_hsc+n_clu);
measured_neg_num_hsc = data_neg_number(:,1+clu_hsc);
errors_neg_num_hsc = data_neg_number(:,1+clu_hsc+n_clu);


measured = [measured_lab_rel(:); measured_cluster_size(:); measured_lab_num_hsc; measured_neg_num_hsc;ones(9,1)*1000];
errors = [errors_lab_rel(:); errors_cluster_size(:); errors_lab_num_hsc; errors_neg_num_hsc;ones(9,1)*500];

%% read M

M = create_differentiation_matrix(cluster_names,n_clu);

%% how many parameters

number_parameters = size(M,1) + (n_clu + 1) *3 +4;

%% optimisation

start = zeros(1,number_parameters);
stop = ones(1,number_parameters) * 4;

stop(1:22) = 200000;
stop(23:44) = 1;
stop(end-1:end) = 500000;
start(21) = 500;
stop(21) = 1500;
start(size(M,1)+(n_clu+1)*2 +3 : end-2) = -4;
stop([36 37]+44) = 0.02;



options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals', 1,'MaxIter',6000);%% change for multistart guess

if exist(name_best, 'file')==2
    
    data2 = dlmread(name_best);
    chibest = data2(end)
    guess = data2(1:number_parameters);
    
else
    
    guess = rand(1,number_parameters) * .01;%% change for multistart guess
    chibest = 1000000000000000000;
    
end



[thetaRes,chisq]=lsqnonlin(@fitFun,guess(:),start,stop,options);



if chisq<chibest
    dlmwrite(name_best,[thetaRes(:);chisq ])
    
    
end


dlmwrite('./output/differentiation_matrix.txt', d, 'delimiter','\t')
dlmwrite('./output/net_proliferation.txt', p, 'delimiter','\t')
dlmwrite('./output/self_renewal.txt', k, 'delimiter','\t')

%% cost function
    function result=fitFun(theta)

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
        
        
        model_lab = deval(sol_lab,time)';
        
        
        model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu + 1) + model_lab(:,n_clu + 2);
        
        
        % solv neg
        
        sol_neg = ode45(@ODE, [time(1), 270], I_neg);
        
        
        model_neg = deval(sol_neg,time)';
        
        model_tip = model_neg(:,n_clu + 1);
        
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
        
        
        model_neg = model_neg ./ repmat(model_neg(:,clu_hsc),1,n_clu);
        model_neg(:,clu_hsc) = [];
        
        
        % ODE system
        
        
        function dxdt = ODE(~,x)
            
            dxdt = d' * x  - k .* x;
            
            xx = x(1)+x(21)+x(22);
            
            dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
            
            
        end
        
        
        
        
        model = [model_lab(:); model_neg(:); model_lab_hsc; model_neg_hsc;model_tip];
        
        
        result = (measured-model)./errors;
        
    end


%% plot

% t_plot = time(1):270;
t_plot = time(1):0.01:270;

%

model_lab = deval(sol_lab, t_plot)';

model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu+1) + model_lab(:,n_clu+2);

model_lab(:,clu_hsc) = model_lab_hsc;

%

model_neg = deval(sol_neg, t_plot)';

model_neg_save = model_neg;


model_tip = model_neg(:,n_clu+1);

model_neg_hsc = model_neg(:,clu_hsc) + model_neg(:,n_clu+1) + model_neg(:,n_clu+2);



model_neg(:,clu_hsc) = model_neg_hsc;



%
model_neg(:,n_clu+1:end) = [];
model_lab(:,n_clu+1:end) = [];


model_lab = model_lab ./ repmat(model_lab(:,clu_hsc),1,n_clu);
model_neg = model_neg ./ repmat(model_neg(:,clu_hsc),1,n_clu);

%% number hsc over time

figure(1001)
clf

subplot(1,3,1)

hold on


errorbar(time,measured_lab_num_hsc,errors_lab_num_hsc, 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)

plot(t_plot,model_lab_hsc,'r','linewidth',1.2)

xlabel('Time (days)')
ylabel('No of labelled cells')

legend('data','model','location','northwest')

subplot(1,3,2)

hold on


errorbar(time,measured_neg_num_hsc,errors_neg_num_hsc, 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)

plot(t_plot,model_neg_hsc,'r','linewidth',1.2)


xlabel('time (d)')
ylabel('Cluster size')

legend('data','model','location','northwest')


subplot(1,3,3)

hold on


% errorbar(time,ones(9,1)*1000,ones(9,1)*500, 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)

plot(t_plot,model_neg_save(:,1),'r','linewidth',1.2)
plot(t_plot,model_neg_save(:,21),'k','linewidth',1.2)
plot(t_plot,model_neg_save(:,22),'c','linewidth',1.2)

legend('0c','0a','0b','location','northwest')


xlabel('Time (days)')
ylabel('Cluster size')

set(gca, 'yscale','log')
ylim([500 inf])

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/real_num_hsc.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width',21.5,'height',4,'Renderer','painters','Lockaxes',0);%


dlmwrite('./output/real_lab.txt',[t_plot(:),model_lab_hsc(:)])
dlmwrite('./output/real_neg.txt',[t_plot(:),model_neg_hsc(:)])
%%

figure(9)
clf

hold on

for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    
    plot(t_plot, model_neg(:,ii), 'r', 'Linewidth', 1.2)
    
    errorbar(time, measured_cluster_size_temp(:,ii), errors_cluster_size_temp(:,ii), 'ob', 'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', 'b', 'MarkerSize', 3)
    
    
    if ii >16 || ii == 16
        xlabel('time (d)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('Rel. cluster size')
    end
    
    title(cluster_names(ii))
    
    
    set(gca,'xscale','log')
end
%
set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/cluster_size.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width',22,'height',25,'Renderer','painters','Lockaxes',0);%

dlmwrite('./output/cluster_size.txt',[t_plot(:),model_neg])

%%
figure(700)
% clf
for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    
    plot(t_plot, model_lab(:,ii), 'r', 'Linewidth', 1.2)
    
    errorbar(time, measured_lab_rel_temp(:,ii), errors_lab_rel_temp(:,ii), 'ob', 'LineWidth', 0.25, 'Capsize', 3, 'MarkerFaceColor', 'b', 'MarkerSize', 3)
    
    %     line([3,300],[1 1],'color','k')
    
    if ii >16 || ii == 16
        xlabel('time (d)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('Rel. label frequency')
    end
    
    title(cluster_names(ii))
    
    
    set(gca,'xscale','log')
end


set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/lab_frequency.eps','FontMode', 'fixed','Fontsize', 6,'color', 'cmyk','width',22,'height',25,'Renderer','painters','Lockaxes',0);%

dlmwrite('./output/lab_freq.txt',[t_plot(:),model_lab])

end
