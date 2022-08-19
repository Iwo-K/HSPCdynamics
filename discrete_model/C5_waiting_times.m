function C5_waiting_times
%% general info

dir = './input/';

name_best = 'output/best.txt';

data2 = dlmread(name_best);

theta = data2(1:104);
theta = theta(:);

cluster_names = [0:12,14,16,20,24,25,26,28];


n_clu = length(cluster_names);
clu_hsc = 0+1;


%% read data


data_lab_rel = dlmread([dir 'input_lab_rel_hsc_pv.txt']);

%% prepare measured vector

time = data_lab_rel(:,1);

%% read M

M = create_differentiation_matrix(cluster_names,n_clu);


%% optimisation

I_neg_temp = theta([1 21 22]);

I_neg = zeros(22,1);

I_neg([1,21,22]) = I_neg_temp/sum(I_neg_temp);
% I_neg([1,21,22]) = I_neg_temp


d = zeros(n_clu+2,n_clu+2);


for index = 1: size(M,1)
    
    d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
    
end



p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);

p = [p(1:20);sum(d(21,:));p(21:end)];% proliferation



r = theta(end-1); % logistic parameter
K = theta(end); % carrying capacity



k = sum(d,2)-p;



% solv neg

sol_neg = ode45(@ODE, [time(1), 270], I_neg);

    function dxdt = ODE(~,x)
        
        
        dxdt = d' * x  - k .* x;
        
        xx = x(1)+x(21)+x(22);
        
        dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));

        
    end

%% plot

t_plot = time(1):0.5:270;

model_neg = deval(sol_neg, t_plot)';

model_neg_hsc = model_neg(:,clu_hsc) + model_neg(:,n_clu+1) + model_neg(:,n_clu+2);

model_neg(:,clu_hsc) = model_neg_hsc;


%
model_neg(:,n_clu+1:end) = [];

indi(1) = 1;
for i = 2:20
    [~,mm]=min(abs(model_neg(:,i)-1));
    
    indi(i) = mm;
end
% [cluster_names(2:end)', t_plot(indi(2:end))']

dlmwrite('./output/waiting_times.txt',[cluster_names(2:end)', t_plot(indi(2:end))'])
%%

figure(9)
clf

hold on

for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    plot(t_plot, model_neg(:,ii), 'r', 'Linewidth', 2)
    
    line([3,300],[1 1])
    line([1,1]*t_plot(indi(ii)),[0.01 1000])
    
    if ii >16 || ii == 16
        xlabel('time (d)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('cl size')
    end
    
    
    
    
    title(cluster_names(ii))
    
    set(gca,'yscale','log')
    
    
    set(gca,'xscale','log')
end

end
