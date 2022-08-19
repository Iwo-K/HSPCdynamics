function C7plot_density

load('./PD_matrixes/data_input_pd_clu_10.mat')
load('./PD_matrixes/results_clu_10.mat')


solbest = 2;


vec = 0:1/298:1;
n_grid = 300;
x = linspace(0,1,n_grid);


%%
figure(222)
clf


subplot(3,3,1)

u1 = ksdensity(D.ind.hist{1,1},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,2},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,3},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,4},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);


hold on


plot(vec,u1(1:299),'k')
plot(vec,u2(1:299),'r')
plot(vec,u3(1:299),'b')
plot(vec,u4(1:299),'m')




total_number(1) = sum(sol{1,solbest}.x(1,:))/299;

plot(vec,sol{1,solbest}.x(1,:)/total_number(1),'k','linewidth',2)


ylabel('density')
xlabel('pseudotime')
title('day 3')
%%



subplot(3,3,2)

u1 = ksdensity(D.ind.hist{1,5},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,6},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,7},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,8},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);
u5 = ksdensity(D.ind.hist{1,9},x,'support',[0,1],'function','pdf');
u5 = u5/trapz(x,u5);
u6 = ksdensity(D.ind.hist{1,10},x,'support',[0,1],'function','pdf');
u6 = u6/trapz(x,u6);

hold on


plot(vec,u1(1:299),'k')
plot(vec,u2(1:299),'r')
plot(vec,u3(1:299),'b')
plot(vec,u4(1:299),'m')
plot(vec,u5(1:299),'g')
plot(vec,u6(1:299),'y')



total_number(2) = sum(sol{1,solbest}.x(2,:))/299;

plot(vec,sol{1,solbest}.x(2,:)/total_number(2),'k','linewidth',2)


ylabel('density')
xlabel('pseudotime')
title('day 7')
%%
subplot(3,3,3)


u1 = ksdensity(D.ind.hist{1,11},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,12},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,13},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,14},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);

hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))


total_number(3) = sum(sol{1,solbest}.x(3,:))/299;

plot(vec,sol{1,solbest}.x(3,:)/total_number(3),'k','linewidth',2)

ylabel('density')
xlabel('pseudotime')
title('day 12')
%%
subplot(3,3,4)


u1 = ksdensity(D.ind.hist{1,15},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,16},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,17},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,18},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);


hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))

total_number(4) = sum(sol{1,solbest}.x(4,:))/299;

plot(vec,sol{1,solbest}.x(4,:)/total_number(4),'k','linewidth',2)


ylabel('density')
xlabel('pseudotime')
title('day 27')
%%

subplot(3,3,5)


u1 = ksdensity(D.ind.hist{1,19},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,20},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,21},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,22},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);


hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))


total_number(5) = sum(sol{1,solbest}.x(5,:))/299;

plot(vec,sol{1,solbest}.x(5,:)/total_number(5),'k','linewidth',2)

ylabel('density')
xlabel('pseudotime')
title('day 49')

%%

subplot(3,3,6)

%
u3 = ksdensity(D.ind.hist{1,23},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,24},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);


hold on


plot(vec,u3(1:299))
plot(vec,u4(1:299))

%
total_number(6) = sum(sol{1,solbest}.x(6,:))/299;

plot(vec,sol{1,solbest}.x(6,:)/total_number(6),'k','linewidth',2)


ylabel('density')
xlabel('pseudotime')
title('day 76')
%%

subplot(3,3,7)


u1 = ksdensity(D.ind.hist{1,25},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,26},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);
u3 = ksdensity(D.ind.hist{1,27},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,28},x,'support',[0,1],'function','pdf');
u4 = u3/trapz(x,u3);

hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))



total_number(7) = sum(sol{1,solbest}.x(7,:))/299;

plot(vec,sol{1,solbest}.x(7,:)/total_number(7),'k','linewidth',2)




ylabel('density')
xlabel('pseudotime')
title('day 112')
%%

subplot(3,3,8)


u1 = ksdensity(D.ind.hist{1,29},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,30},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);


u3 = ksdensity(D.ind.hist{1,31},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,32},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);

hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))


total_number(8) = sum(sol{1,solbest}.x(8,:))/299;

plot(vec,sol{1,solbest}.x(8,:)/total_number(8),'k','linewidth',2)



ylabel('density')
xlabel('pseudotime')
title('day 161')
%%
subplot(3,3,9)


u1 = ksdensity(D.ind.hist{1,33},x,'support',[0,1],'function','pdf');
u1 = u1/trapz(x,u1);
u2 = ksdensity(D.ind.hist{1,34},x,'support',[0,1],'function','pdf');
u2 = u2/trapz(x,u2);


u3 = ksdensity(D.ind.hist{1,35},x,'support',[0,1],'function','pdf');
u3 = u3/trapz(x,u3);
u4 = ksdensity(D.ind.hist{1,36},x,'support',[0,1],'function','pdf');
u4 = u4/trapz(x,u4);

hold on


plot(vec,u1(1:299))
plot(vec,u2(1:299))
plot(vec,u3(1:299))
plot(vec,u4(1:299))

total_number(9) = sum(sol{1,solbest}.x(9,:))/299;
% 
plot(vec,sol{1,solbest}.x(9,:)/total_number(9),'k','linewidth',2)


ylabel('density')
xlabel('pseudotime')
title('day 269')
%%

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/density_clu_10.eps','FontMode', 'fixed','Fontsize',13,'color', 'cmyk','width',30,'height',20,'Renderer','painters','Lockaxes',0);%

dlmwrite('./tables/total_number_clu_10.txt',total_number)

end
