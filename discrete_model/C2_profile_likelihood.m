function C2_profile_likelihood

name_best = './output/best.txt';
data2 = dlmread(name_best);
dir = './input/';

%% general info

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

%% code

up = 1;
down = 1;

steps = 1000;
countmax = 1;
threshold = 3.8;

small_steps = 0.00001;

options = optimset('display','iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',1500,'MaxIter',6000);%



%%
chibest = data2(end);

best = data2(1:number_parameters);


chimin = chibest;

best_bu = best;

%%

start_base = zeros(1,number_parameters);
stop_base = ones(1,number_parameters) * 4;


stop_base(1:22) = 200000;
stop_base(23:44) = 1;
stop_base(end-1:end) = 500000;
start_base(21) = 500;
stop_base(21) = 1500;
start_base(size(M,1)+(n_clu+1)*2 +3 : end-2) = -4;
stop_base([36 37]+44) = 0.02;


for lookAt = 1:number_parameters
    
    start = start_base;
    start(lookAt)=[];
    
    stop = stop_base;
    stop(lookAt)=[];

    fName1 = ['./output/profile_likelihood/r_',num2str(lookAt),'tot.txt'];
    
    fName2 = ['./output/profile_likelihood/r_',num2str(lookAt),'.txt'];
    
    
    boolFixed=zeros(1,number_parameters);
    boolFixed(lookAt)=1;
    boolFixed=logical(boolFixed);
    
    thetaStart=best;
    
    
    if down
        
        if exist(fName1,'file')==2
            
            
            data3 = dlmread(fName1);
            
            data3(data3(:,end)>chibest+3.8,:) = [];
            
            column = data3(:,lookAt);
            
            
            
            [~, index]=min(column);
            
            best=data3(index,1:number_parameters);
            
            
        end
        

        bestFitParameter = best(lookAt);
        
        step_down = abs(bestFitParameter)/steps;
        
        count = 0;
        
        if  abs(bestFitParameter)<0.001
            
            
            step_down = small_steps;
        end
        
        
        
        bestFitParameter = bestFitParameter-step_down;
        
        while count<countmax && bestFitParameter>=start_base(lookAt)
            
            
            
            thetaStart = best';
            
            
            thetaStart(lookAt)=bestFitParameter;
            

            [reducedResult,chisq]=lsqnonlin(@fitFun,thetaStart(~boolFixed),start,stop,options);
            bestFitParameter= bestFitParameter-step_down;
            theta=zeros(1,number_parameters);
            theta(boolFixed)=thetaStart(boolFixed);

            
            theta(~boolFixed) = reducedResult;
            
            
            
            if chisq>chibest+threshold
                count=count+1;
                
                
            else
                step_down = step_down *1.01;
                count=0;
                best=theta;
                dlmwrite(fName2,[theta(lookAt),chisq],'-append');
                dlmwrite(fName1,[theta,chisq],'-append');
            end
            
            if chisq<chimin
                chimin=chisq;
                
                
            end
            
        end
        
    end
    
    
    if up
        
        best=best_bu;
        
        start=start_base;
        start(lookAt) = [];
        
        stop=stop_base;
        stop(lookAt)=[];
        
        if exist(fName1,'file')==2
            
            
            data3=dlmread(fName1);
            data3(data3(:,end)>chibest+3.8,:) = [];
            
            
            column=data3(:,lookAt);
            
            [~, index]=max(column);
            
            best=data3(index,1:number_parameters);
            
            
        end

        
        bestFitParameter=   best(lookAt);
        
        
        
        step_up = abs(bestFitParameter)/steps;
        if  abs(bestFitParameter)<0.001
            
            
            step_up=small_steps;
        end
        
        
        
        count=0;
        
        bestFitParameter= bestFitParameter+step_up;
        while count<countmax && bestFitParameter<=stop_base(lookAt)
            
            
            
            thetaStart=best;
            

            thetaStart(lookAt)=bestFitParameter;
            
            [reducedResult,chisq]=lsqnonlin(@fitFun,thetaStart(~boolFixed),start,stop,options);
            bestFitParameter= bestFitParameter+step_up;
            
            
            theta=zeros(1,number_parameters);
            theta(boolFixed)=thetaStart(boolFixed);
            theta(~boolFixed)=reducedResult;
            
            
            if chisq>chibest+threshold
                count=count+1;
                
                
            else
                
                step_up = step_up *1.01;
                
                dlmwrite(fName2,[theta(lookAt),chisq],'-append');
                dlmwrite(fName1,[theta,chisq],'-append');
                
                count=0;
                best=theta;
            end
            
            if chisq<chimin
                chimin=chisq;
                
                
            end
            
        end
        
        
        
    end
    
    
end

    function result=fitFun(reducedtheta)
        
        
        theta=zeros(number_parameters,1);
        theta(boolFixed)=thetaStart(boolFixed);
        theta(~boolFixed)=reducedtheta;
        
        
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
        
        sol_lab = ode113(@ODE, [time(1), 270], l0);
        
        
        model_lab = deval(sol_lab,time)';
        
        
        
        model_lab_hsc = model_lab(:,clu_hsc) + model_lab(:,n_clu + 1) + model_lab(:,n_clu + 2);
        
        
        % solv neg
        
        sol_neg = ode113(@ODE, [time(1), 270], I_neg);
        
        
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

end
