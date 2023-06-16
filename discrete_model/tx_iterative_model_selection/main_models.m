function main_models(vec_para,store_changes)


input_dir = './input_data/';

time = 3:0.1:7;

%% general info
vec_para = unique(vec_para);

vec_para = sort(vec_para);


vec_dif = sort(vec_para(vec_para<=37))


ldif = length(vec_dif);
vec_prol = sort(vec_para((vec_para>37) & (vec_para<=58))-37)

lprol = length(vec_prol);

vec_log = sort(vec_para(vec_para>58)-58)


%%

num_dec_dif = sum(2.^(37-vec_dif));


num_hexa_dif = dec2hex(num_dec_dif);


num_dec_prol = sum(2.^(21-vec_prol));


num_hexa_prol = dec2hex(num_dec_prol);


num_dec_log = sum(2.^(2-vec_log));


num_hexa_log = dec2hex(num_dec_log);


%%


name_best = ['./bests/best_' num_hexa_dif '_' num_hexa_prol '_' num_hexa_log '.txt']

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;


%% read data

data = readtable([input_dir 'tx_cellnos.csv']);
data_m = [data.x1Day, data.x3Day, data.x5Day, data.x7Day]';
data_norm = data_m ./ repmat(data_m(:,clu_hsc),1,n_clu);% rel_increase = data_calr.mean_ratio;
% err_data_norm = dlmread([input_dir './input_tx_boot_err_prop.txt']);
err_data_norm = repmat(0.1*max(data_norm),4,1);

err_data_norm(err_data_norm == 0) = 0.1;

% size(data_norm)
% size(err_data_norm)
% return



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



%% prepare measured vector

%% read M

M = create_differentiation_matrix(cluster_names,n_clu);

%% how many parameters

number_parameters = length(vec_para)
%% optimisation



options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals', 100000,'MaxIter',6000);%% change for multistart guess

if exist(name_best, 'file')==2

    data2 = dlmread(name_best);
    chibest = data2(end)

    guess = data2(1:number_parameters);

else

    guess = rand(1,number_parameters) * .001;%% change for multistart guess
    chibest = 1000000000000000000;

end

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




start = -theta_base(44+vec_para);
stop = 4-theta_base(44+vec_para);


start(ldif+1:ldif+lprol) = -4-theta_base(44+37+vec_prol);

vv = [0.4,200000];

stop(ldif+lprol+1:end) = ones(1,number_parameters-lprol-ldif).*vv(vec_log);
start(ldif+lprol+1:end) = -ones(1,number_parameters-lprol-ldif).*vv(vec_log);



try

    chisq = 0;
    delta = chibest - chisq;


    while delta > 1

        [thetaRes,chisq]=lsqnonlin(@fitFun,guess(:),start,stop,options);



        if chisq<chibest
            dlmwrite(name_best,[thetaRes(:);chisq ])
        end


        delta = chibest - chisq;
        chibest = chisq;
        guess = thetaRes;

    end
    %
catch
    'no best'
end





%%
    function dxdt = ODE(~,x)

        dxdt = d' * x  - k .* x;



        xx = x(1)+x(21)+x(22);

        dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
        %

    end

%% cost function
    function result=fitFun(theta)
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



        % solv neg
        try


            sol = ode45(@ODE, [time(1), 8], c3);



            model3_best = deval(sol,[3 5 7])';
            % %Re-normalise to cluster 0 (sum of cluster 0,30,40)

            hsc = model3_best(:,clu_hsc) + model3_best(:,n_clu+1) + model3_best(:,n_clu+2);
            model3_best(:,clu_hsc) = hsc;


            model3_best(:,n_clu+1:end) = [];
            model3_best = model3_best ./ repmat(model3_best(:,clu_hsc),1,n_clu);



            %
        catch
            'error'

        end


        % ODE system


        function dxdt = ODE(~,x)

            dxdt = d' * x  - k .* x;



            xx = x(1)+x(21)+x(22);

            dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
            %

        end


        result = (data_norm-model3_best)./err_data_norm;
    end


%%


end




