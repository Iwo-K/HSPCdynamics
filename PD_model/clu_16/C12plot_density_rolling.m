function C12plot_density_rolling

data_roll =  load('./PD_matrixes/results_rolling_clu_16.mat');



solbest = 3;


vec = 0:1/298:1;


%%

k = 0;

figure(1)
clf

hold on


for n = [ 1 20 40 51]
    
    
    hold on
    
k = k+1;
    
    
    l(k) = plot(vec,data_roll.sol{1,solbest}.x(n,:),'linewidth',2);
    
    
    ylabel('density')
    xlabel('pseudotime')
    
    
end



nn = n

yy = data_roll.sol{1,solbest}.x(nn,:);

dpt = 0.6;

AUC = trapz(vec(vec >= dpt),yy(vec >= dpt))



vvx = vec(vec>=dpt);

vvyt = data_roll.sol{1,solbest}.x(nn,:);
vvy = vvyt(vec>=dpt);

X = [vvx'; flipud(vvx')]';


Y = [vvy'; flipud(vvx')]';




l(end+1) = fill(X, Y, 'b', 'facecolor', [1 1 1] * 0.8, 'edgecolor', [1 1 1] * 0.8);


legend(l,'1 day','20 days','40 days','51 days','1 cell arrived')




set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/waiting_times_rolling_clu_16.eps','FontMode', 'fixed','Fontsize',13,'color', 'cmyk','width',12,'height',10,'Renderer','painters','Lockaxes',0);%

end

