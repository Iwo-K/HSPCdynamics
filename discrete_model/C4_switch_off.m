function C4_switch_off(vec_clu)
%% general info

cluster_names = [0:12,14,16,20,24,25,26,28,30,40];

n_clu = length(cluster_names)-2;
clu_hsc = 0+1;


name_best = './output/best.txt';

[~,vec] = ismember(vec_clu,  cluster_names);


%% read M

M = create_differentiation_matrix(cluster_names,n_clu);

%% how many parameters

number_parameters = size(M,1) + (n_clu + 1) *3 +4;

%% generate curves



data2 = dlmread(name_best);

theta = data2(1:number_parameters);
theta = theta(:);


%

I_neg = theta(1: n_clu+2);



d = zeros(n_clu+2,n_clu+2);% differentiation


for index = 1: size(M,1)
    
    d(M(index,1),M(index,2)) = theta(index+(n_clu+1)*2 +2);
    
end



p = theta(size(M,1)+(n_clu+1)*2 +3 : end-2);

p = [p(1:20);sum(d(21,:));p(21:end)];% proliferation



r = theta(end-1); % logistic parameter
K = theta(end); % carrying capacity



k = sum(d,2)-p;





% solv tot

sol_tot_1 = ode45(@ODE, [1, 270], I_neg);


% solv tot off
I_neg(vec_clu) = 0;

sol_tot_0 = ode45(@ODE, [1, 270], I_neg);




% ODE system

    function dxdt = ODE(~,x)
        
        
        dxdt = d' * x  - k .* x;
        
        xx = x(1)+x(21)+x(22);
        
        dxdt(21) = r*xx*(1-xx/K) - (dxdt(1)+dxdt(22));
        
    end






%% plot

t_plot = 1:270;





model_tot1 = deval(sol_tot_1, t_plot)';


model_tot_hsc1 = model_tot1(:,clu_hsc) + model_tot1(:,n_clu+1) + model_tot1(:,n_clu+2);

model_tot1(:,clu_hsc) = model_tot_hsc1;



model_tot0 = deval(sol_tot_0, t_plot)';


model_tot_hsc0 = model_tot0(:,clu_hsc) + model_tot0(:,n_clu+1) + model_tot0(:,n_clu+2);

model_tot0(:,clu_hsc) = model_tot_hsc0;

model_R = model_tot0./model_tot1;
%% cluster size

figure(9)


hold on

for ii = 1:n_clu
    
    subplot(5,4,ii)
    hold on
    
    plot(t_plot, model_R(:,ii), 'r', 'Linewidth', 2)
    
    
    
    if ii >16 || ii == 16
        xlabel('time (d)')
        
    end
    if ii == 1 ||  ii == 5 ||  ii == 9 ||  ii == 13 ||  ii == 17
        ylabel('relative cluster size')
    end
    
    title(cluster_names(ii))
    
    
    set(gca,'xscale','log')
end


dlmwrite(['./output/waiting_times_deterministic' num2str(vec_clu) '.txt'],[cluster_names', model_R']);


set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,['./figures/rel_clu_size_after_switching_off',num2str(vec_clu),'.eps'],'FontMode', 'fixed','Fontsize',13,'color', 'cmyk','width',50,'height',50,'Renderer','painters','Lockaxes',0);%


end