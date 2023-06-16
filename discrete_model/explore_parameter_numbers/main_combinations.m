function main_loop(vec)
%% general info

dir = './input/';

name_best_base = [dir  'best_base.txt'];


name_best = ['./output/combinations/best_' num2str(dec2hex(bin2dec(num2str(vec)))) '.txt']

vec_i = find(vec);


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


if exist(name_best, 'file')==2


    'already there'
    return
else
    for ee = 10:-1:1
        ee
        measured = [measured_lab_rel(:); measured_cluster_size(:); measured_lab_num_hsc; measured_neg_
            num_hsc;ones(9,1)*1000];
        errors = [errors_lab_rel(:); errors_cluster_size(:); errors_lab_num_hsc/ee; errors_neg_num_hsc/ee;ones(9,1)*500/ee];

        %% read M

        M = create_differentiation_matrix(cluster_names,n_clu);
        M = M(vec_i,:);


        %% how many parameters

        number_parameters_b = 104;
        number_parameters = number_parameters_b + length(vec_i) - 37

        %% optimisation

        start_base = zeros(1,number_parameters_b);
        stop_base = ones(1,number_parameters_b) * 4;

        stop_base(1:22) = 200000;
        stop_base(23:44) = 1;
        stop_base(end-1:end) = 500000;
        start_base(21) = 500;
        stop_base(21) = 1500;
        start_base(37+(n_clu+1)*2 +3 : end-2) = -4;
        stop_base([36 37]+44) = 0.02;

        start = start_base([1:44,44+vec_i,number_parameters_b-1-21:number_parameters_b]);
        stop = stop_base([1:44,44+vec_i,number_parameters_b-1-21:number_parameters_b]);

        %

        options=optimset('Display','Iter','TolX',1e-8,'TolFun',1e-8,'MaxFunEvals', round(10000/ee),'MaxIter',6000);%% change for multistart guess

        if exist(name_best, 'file')==2

            data2 = dlmread(name_best);
            chibest = data2(end)
            guess = data2(1:number_parameters);

        else

            data2 = dlmread(name_best_base);

            best_base = data2(1:number_parameters_b);

            guess = best_base([1:44,44+vec_i,number_parameters_b-1-21:number_parameters_b]);
            % sum(guess(45:44+37)<0.00001)

            %     guess(guess>=0.1) = 0.1;
            guess = rand(number_parameters,1)*0.0001;
            chibest = 10^100;

        end

        value_continue = 1;

        while value_continue

            try
                [thetaRes,chisq]=lsqnonlin(@fitFun,guess(:),start,stop,options);
            catch
                'did not converge'
                return
            end


            if chisq<chibest-1
                dlmwrite(name_best,[thetaRes(:);chisq ])

                data2 = dlmread(name_best);
                chibest = data2(end)
                guess = data2(1:number_parameters);


            else
                value_continue = 0;

            end

        end

    end
end

%% cost function
    function result=fitFun(theta)

        I_neg = theta(1: n_clu+2);


        l0 = theta(n_clu+3:(n_clu+1)*2 +2).*I_neg;

        d = zeros(n_clu+2,n_clu+2);% differentiation


        for index = 1: size(M,1)

            d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);

        end



        p = theta(44+length(vec_i)+1 : end-2);

        p = [p(1:20);sum(d(21,:));p(21:end)];% proliferation


        r = theta(end-1); % logistic parameter
        K = theta(end); % carrying capacity


        k = sum(d,2)-p;


        % solv lab

        try
            sol_neg = ode45(@ODE, [time(1), 270], I_neg);
        catch
            'problem using ode45'
            return
        end

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
