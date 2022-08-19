function C12plot_density_rolling

data_roll =  load('./PD_matrixes/results_rolling_clu_10.mat')



solbest = 2;


vec = 0:1/298:1;


%%



figure(1)
clf

hold on

kk = 0;
for n = [ 1 10 20 36]
    
    
    hold on
    
    
    kk = kk +1
    
    l(kk) = plot(vec,data_roll.sol{1,solbest}.x(n,:),'linewidth',2);
    
    ylabel('density')
    xlabel('pseudotime')
    
    
end

nn = n;
yy = data_roll.sol{1,solbest}.x(nn,:);

AUC = trapz(vec(vec >= 0.6),yy(vec >= 0.6))



vvx = vec(vec>=0.6);

vvyt = data_roll.sol{1,solbest}.x(nn,:);
vvy = vvyt(vec>=0.6);

X = [vvx'; flipud(vvx')]';


Y = [vvy'; flipud(vvx')]';




l(kk+1) = fill(X, Y, 'b', 'facecolor', [1 1 1] * 0.8, 'edgecolor', [1 1 1] * 0.8)


legend(l,'1 day','10 days','20 days','36 days','1 cell arrived');

%%

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/waiting_times_rolling_clu_10.eps','FontMode', 'fixed','Fontsize',13,'color', 'cmyk','width',12,'height',10,'Renderer','painters','Lockaxes',0);%


end

