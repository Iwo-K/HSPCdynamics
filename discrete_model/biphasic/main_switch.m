function main_switch
%% general info

dir = '../input/';

name_best = './output/best_switch.txt';

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;

sp = 4; %switching point
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

number_parameters = size(M,1) *2 + (n_clu + 2) *2 + (n_clu+1)*2 +4;

%% optimisation

start = zeros(1,number_parameters);
stop = ones(1,number_parameters) * 4;

stop(1:22) = 200000;
stop(23:44) = 1;
stop(103:104) = 5000000;
start(21) = 500;
stop(21) = 1500;
start(82:102) = -4;
stop([36 37]+44) = 0.02;


stop(163:164) = 5000000;
start(142:162) = -4;
stop([36 37]+44+60) = 0.02;


options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals', 5,'MaxIter',6000);%% change for multistart guess

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

%% cost function
    function result=fitFun(theta)

        
        I_neg = theta(1:22);
        
        
        l0 = theta(23:44) .* I_neg;
        
        d = zeros(n_clu+2,n_clu+2);% differentiation
        
        
        for index = 1: size(M,1)
            
            d(M(index,1),M(index,2)) = theta(index + 44);
            
        end
        
        
        
        p = theta(82:102);
        
        p = [p(1:20);sum(d(21,:));p(21)];% proliferation
        
        
        
        r = theta(103); % logistic parameter
        K = theta(104); % carrying capacity
        
        
        k = sum(d,2)-p;
        
        
        % solv lab
        
        sol_lab1 = ode45(@ODE, [time(1), time(sp)], l0);
        
        
        model_lab1 = deval(sol_lab1,time(1:sp))';
        
        
        
        % solv neg
        
        sol_neg1 = ode45(@ODE, [time(1), time(sp)], I_neg);
        
   
        model_neg1 = deval(sol_neg1,time(1:sp))';
        
        
        % 2

        d = zeros(n_clu+2,n_clu+2);% differentiation
        
        
        for index = 1: size(M,1)
            
            d(M(index,1),M(index,2)) = theta(index+104);
            
        end
        
        
        p = theta(142:162);
        
        p = [p(1:20);sum(d(21,:));p(21:end)];% proliferation
        
        
        
        r = theta(163); % logistic parameter
        K = theta(164); % carrying capacity
        
        
        k = sum(d,2)-p;
        
        %

        
        l02 = model_lab1(end,:);
        
        
        sol_lab2 = ode45(@ODE, [time(sp), 270], l02);
        
        
        model_lab2 = deval(sol_lab2,time(sp+1:end))';
        
        
        
        model_lab = [model_lab1;model_lab2];
        
        
        model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu + 1) + model_lab(:,n_clu + 2);
        
        
    
        
        
        I_neg2 = model_neg1(end,:);
        
        
        sol_neg2 = ode45(@ODE, [time(sp), 270], I_neg2);
        
        
        model_neg2 = deval(sol_neg2,time(sp+1:end))';
        
        
        
        model_neg = [model_neg1;model_neg2];
        
        
        
        
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

t_plot1 = time(1):0.01:time(sp);
t_plot2 = time(sp):0.01:270;
t_plot = [t_plot1,t_plot2];
t_plot_one_phase = time(1):0.01:270;

%

model_lab1 = deval(sol_lab1, t_plot1)';
model_lab2 = deval(sol_lab2, t_plot2)';

model_lab = [model_lab1;model_lab2];

model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu+1) + model_lab(:,n_clu+2);

model_lab(:,clu_hsc) = model_lab_hsc;

%

model_neg1 = deval(sol_neg1, t_plot1)';
model_neg2 = deval(sol_neg2, t_plot2)';

model_neg = [model_neg1;model_neg2];


model_neg_save = model_neg;

model_tip = model_neg(:,n_clu+1);

model_neg_hsc = model_neg(:,clu_hsc) + model_neg(:,n_clu+1) + model_neg(:,n_clu+2);





model_neg(:,clu_hsc) = model_neg_hsc;



%
model_neg(:,n_clu+1:end) = [];
model_lab(:,n_clu+1:end) = [];


model_lab = model_lab ./ repmat(model_lab(:,clu_hsc),1,n_clu);
model_neg = model_neg ./ repmat(model_neg(:,clu_hsc),1,n_clu);


%%

data_one_phase = dlmread('../output/cluster_size.txt');
model_neg_one_phase = data_one_phase(:,2:end);


figure(9)
clf

hold on
col = [0.55, 0.55, 0.55]
for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    
    errorbar(time, measured_cluster_size_temp(:,ii), errors_cluster_size_temp(:,ii), 'ob', 'Color', col,  'LineWidth', 0.25, 'Capsize', 3,  'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 3)

    l1 = plot(t_plot, model_neg(:,ii), 'b', 'Linewidth', 1.2)
    l2 = plot(t_plot_one_phase, model_neg_one_phase(:,ii), 'r', 'Linewidth', 1.2)

    
    
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
exportfig(gcf,'./figures/biphasic_fit.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width', 32,'height',25,'Renderer','painters','Lockaxes',0);%




end
