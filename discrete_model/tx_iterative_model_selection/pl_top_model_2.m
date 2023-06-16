function pl_top_model_2

vec_para = [3,4,6,7,8,9,10,11,12,13,16,18,19,20,23,24,25,26,27,31,32,33,39,40,46,50,54];

input_dir = './input_data/';

time = 3:0.1:7;
%% general info

vec_para = unique(vec_para);


vec_para = sort(vec_para);



vec_dif = sort(vec_para(vec_para<=37));



ldif = length(vec_dif);
vec_prol = sort(vec_para((vec_para>37) & (vec_para<=58))-37);

lprol = length(vec_prol);

vec_log = sort(vec_para(vec_para>58)-58);

%%

num_dec_dif = sum(2.^(37-vec_dif));


num_hexa_dif = dec2hex(num_dec_dif);


num_dec_prol = sum(2.^(21-vec_prol));


num_hexa_prol = dec2hex(num_dec_prol);


num_dec_log = sum(2.^(2-vec_log));


num_hexa_log = dec2hex(num_dec_log);

%%
name_best = ['./bests/best_' num_hexa_dif '_' num_hexa_prol '_' num_hexa_log '.txt']
data2 = dlmread(name_best);

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;

%% read data

data = readtable([input_dir 'tx_cellnos.csv']);
data_m = [data.x1Day, data.x3Day, data.x5Day, data.x7Day]';
data_norm = data_m ./ repmat(data_m(:,clu_hsc),1,n_clu);% rel_increase = data_calr.mean_ratio;
err_data_norm = repmat(0.1*max(data_norm),4,1);
err_data_norm(err_data_norm == 0) = 0.1;

data_norm(1,:) = [];
err_data_norm(1,:) = [];

% Loading cluster sizes
best = dlmread([input_dir 'best_base.txt']);
best(end) = [];
clu_size = best(1:22);
hsc_totsize = sum(clu_size([clu_hsc n_clu+1 n_clu+2]));
cl0_relsize = clu_size(clu_hsc) / hsc_totsize;
cl30_relsize = clu_size(n_clu + 1) / hsc_totsize;
cl40_relsize = clu_size(n_clu + 2) / hsc_totsize;
%% Day3 as starting point

c3 = data.x3Day;
c3 = c3/c3(1);
% Splitting cluster 0 into: 0, 30, 40 with in the proportions from the Hoxb5
% model
c3(clu_hsc) = cl0_relsize;
c3 = [c3; cl30_relsize; cl40_relsize];


%% read M

M = create_differentiation_matrix(cluster_names,n_clu);

%% how many parameters

number_parameters = length(vec_para)

%% code

up = 1;
down = 1;

steps = 8103;
countmax = 1;
threshold = 3.8;

small_steps = 0.00001;

options = optimset('display','iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',25000,'MaxIter',6000);%



%%
chibest = data2(end);

best = data2(1:number_parameters);


chimin = chibest;

best_bu = best;

%%

theta_base = dlmread('./input_data/best_base.txt');


theta_dif = theta_base(45:44+37);

% length(theta_dif)

theta_mix = theta_base;
theta_mix(end) = [];
theta_base(end) = [];

r = theta_mix(end-1); % logistic parameter
K = theta_mix(end); % carrying capacity


r_base = r;
K_base = K;


%%

start_base = -theta_base(44+vec_para);
stop_base = 4-theta_base(44+vec_para);


start_base(ldif+1:ldif+lprol) = -4-theta_base(44+37+vec_prol);

vv = [0.4,200000];

stop_base(ldif+lprol+1:end) = ones(1,number_parameters-lprol-ldif).*vv(vec_log);
start_base(ldif+lprol+1:end) = -ones(1,number_parameters-lprol-ldif).*vv(vec_log);



for lookAt = 1:27
    
    start = start_base;
    start(lookAt)=[];
    
    stop = stop_base;
    stop(lookAt)=[];

    fName1 = ['./pl_tm_2/r_',num2str(lookAt),'tot.txt'];
    
    fName2 = ['./pl_tm_2/r_',num2str(lookAt),'.txt'];
    
    
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
        
        
        theta = theta(:);



        theta_dif(vec_dif) = theta(1:ldif) + theta_base(44+ vec_dif);



        d = zeros(n_clu+2,n_clu+2);% differentiation


        for index = 1: size(M,1)
            d(M(index,1),M(index,2)) = theta_dif(index);

        end

        p = theta_base(size(M,1)+44+1 : end-2);


        p(vec_prol) = theta(ldif+1:ldif+lprol) + p(vec_prol);

        p = [p(1:20);sum(d(21,:));p(21)];




        k = sum(d,2)-p;

        if ismember(59,vec_para) && ismember(60,vec_para)

            r = theta(ldif+lprol+1) + r_base;
            K = theta(ldif+lprol+2)+ K_base;

        elseif ismember(59,vec_para)
            r = theta(ldif+lprol+1)+ r_base;

        elseif ismember(60,vec_para)
            K = theta(ldif+lprol+1) + K_base;
        end




        sol = ode45(@ODE, [time(1), 8], c3);



        model3_best = deval(sol,[3 5 7])';
        % %Re-normalise to cluster 0 (sum of cluster 0,30,40)

        hsc = model3_best(:,clu_hsc) + model3_best(:,n_clu+1) + model3_best(:,n_clu+2);
        model3_best(:,clu_hsc) = hsc;
        model3_best(:,n_clu+1:end) = [];
        model3_best = model3_best ./ repmat(model3_best(:,clu_hsc),1,n_clu);





        % ODE system


        function dxdt = ODE(~,x)

            dxdt = d' * x  - k .* x;



            xx = x(1)+x(21)+x(22);

            dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
            %

        end


        result = (data_norm-model3_best)./err_data_norm;
        
    end

end