function C1_bootstrap(n)
%% general info
rng(n*now)

name_best = '../output/best.txt';

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;

%% read data

data_lab_rel = create_input_bootstrap_lab_rel;
data_lab_number = create_input_bootstrap_lab_flow_nodivide;

data_cluster_size = create_input_bootstrap_cluster_size_relative_to_hsc;
data_neg_number = create_input_bootstrap_neg_flow_nodivide;


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


qq = 1;
measured = [measured_lab_rel(:); measured_cluster_size(:); measured_lab_num_hsc; measured_neg_num_hsc;ones(9,1)*1000];
errors = [errors_lab_rel(:); errors_cluster_size(:); errors_lab_num_hsc/qq; errors_neg_num_hsc/qq;ones(9,1)*500/qq];
%% read M


M = create_differentiation_matrix(cluster_names,n_clu);

%% how many parameters

number_parameters = size(M,1) + (n_clu + 1) *3 +4

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


options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals', 300000,'MaxIter',60000);


data2=dlmread(name_best);

guess=data2(1:number_parameters);


[theta,chisq,~,flag]=lsqnonlin(@fitFun,guess(:),start,stop,options);

dlmwrite('bootstrap_simulations.txt',[theta;chisq;flag]','-append');

%% cost function
    function result=fitFun(theta)
        
        
        I_neg = theta(1: n_clu+2);
        
        
        l0 = theta(n_clu+3:(n_clu+1)*2 +2).*I_neg;
        
        d = zeros(n_clu+2,n_clu+2);
        
        
        for index = 1: size(M,1)
            
            d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
            
        end
        
        
        p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);
        
        p = [p(1:20);sum(d(21,:));p(21:end)];
        
        
        r = theta(end-1);
        K = theta(end);
        

        k = sum(d,2)-p;
        
        
        % solv lab
        
        sol_lab = ode45(@ODE, [time(1), 270], l0);
        
        %     [X,Y] = ode45(@ODE, [time(1), 270], l0);
        
        
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
        
        function dxdt = ODE(~,x)

            dxdt = d' * x  - k .* x;
            
            xx = x(1)+x(21)+x(22);
            
            dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
            
        end
        
        
        model = [model_lab(:); model_neg(:); model_lab_hsc; model_neg_hsc;model_tip];
        
        
        result = (measured-model)./errors;
        
        
    end


end
