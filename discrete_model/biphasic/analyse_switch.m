function analyse_switch

data = dlmread('./output/profile_likelihood/bounds_switch.txt');


data(1:44,:) = [];


%% differentiation

label={'0b-4',...
    '0b-5',...
    '0b-8',...
    '0c-8',...
    '0b-0c',...
    '1-9',...
    '1-12',...
    '2-3',...
    '2-4',...
    '2-5',...
    '2-6',...
    '2-12',...
    '2-16',...
    '3-10',...
    '3-25',...
    '4-2',...
    '4-5',...
    '4-8',...
    '4-12',...
    '5-2',...
    '5-4',...
    '5-16',...
    '6-24',...
    '8-1',...
    '8-4',...
    '8-7',...
    '8-12',...
    '9-11',...
    '10-25',...
    '11-20',...
    '12-25',...
    '12-26',...
    '24-14',...
    '16-14',...
    '16-28',...
    '0a-0b',...
    '0a-0c'};




figure(9)
clf

hold on

for i = 1:37
    
    best = data(i,2);
    xneg = best - data(i,4);
    xpos = data(i,3) - best;
    
    l1 = plot(best,i,'ob','markerfacecolor', 'b', 'MarkerSize', 2);
    
    errorbar(best,i,xneg,xpos,'horizontal', 'color','b', 'LineWidth', 0.15, 'Capsize', 1,  'MarkerFaceColor', 'b', 'MarkerSize', 0.25)
    
    best = data(i+60,2);
    xneg = best - data(i+60,4);
    xpos = data(i+60,3) - best;
    
    l2 = plot(best,i,'or','markerfacecolor','r', 'MarkerSize', 2);
    
    errorbar(best,i,xneg,xpos,'horizontal','color','r', 'LineWidth', 0.15, 'Capsize', 1,  'MarkerFaceColor', 'r', 'MarkerSize', 0.25)
end


set(gca,'ytick',[])
set(gca,'yTick',1:37,'Box','off')
set(gca,'yTickLabel',label)

xlabel('Differentiation day^{-1}')
ylabel('Cluster-cluster transition')


axis([-0.1 4.1 0.3 37.2])

legend([l1,l2],'Phase I','Phase II')

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/diff_rates_biphasic.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width', 8,'height',11,'Renderer','painters','Lockaxes',0);%



%% proliferation

label={'0b',...
    '1',...
    '2',...
    '3',...
    '4',...
    '5',...
    '6',...
    '7',...
    '8',...
    '9',...
    '10',...
    '11',...
    '12',...
    '14',...
    '16',...
    '20',...
    '24',...
    '25',...
    '26',...
    '28',...
    '0c',...
    };


figure(1)
clf

hold on

for i = 38:58
    
    best = data(i,2);
    xneg = best - data(i,4);
    xpos = data(i,3) - best;
    
    l1 = plot(best,i-37,'ob','markerfacecolor','b', 'MarkerSize', 2);
%     s1 = scatter(best,i-37, 20,'blue', 'filled');
%     s1 = scatter(best,i-37, 20, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    
    e1 = errorbar(best,i-37,xneg,xpos,'horizontal','color','b', 'LineWidth', 0.15, 'Capsize', 1,  'MarkerFaceColor', 'b', 'MarkerSize', 0.25);
%     set([e1.Bar, e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*0.5]);

    best = data(i+60,2);
    xneg = best - data(i+60,4);
    xpos = data(i+60,3) - best;
    
    l2 = plot(best,i-37,'or','markerfacecolor','r', 'MarkerSize', 2);
%     s2 = scatter(best,i-37, 20,'red', 'filled');
%     s2 = scatter(best,i-37, 20, 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

    e2 = errorbar(best,i-37,xneg,xpos,'horizontal','color','r', 'LineWidth', 0.15, 'Capsize', 1,  'MarkerFaceColor', 'r', 'MarkerSize', 0.25);
%     set([e2.Bar, e2.Line], 'ColorType', 'truecoloralpha', 'ColorData',
%     [e2.Line.ColorData(1:3); 255*0.5]); 
end


set(gca,'ytick',[])
set(gca,'yTick',1:21,'Box','off')
set(gca,'yTickLabel',label)
axis([-0.1 4.1 0.3 37.2])
legend([l1,l2],'Phase I','Phase II')
axis([-4.1 4.1 0.3 21.2])

xlabel('Net proliferation day^{-1}')
ylabel('Cluster')

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/netprolif_rates_biphasic.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width', 8,'height',8,'Renderer','painters','Lockaxes',0);%


return

%% the other 2 rates ???

figure(2)
clf

hold on

for i = 59:60
    
    best = data(i,2);
    xneg = best - data(i,4);
    xpos = data(i,3) - best;
    
    plot(best,i-58,'ob','markerfacecolor','b')
    
    errorbar(best,i-58,xneg,xpos,'horizontal','color','b')
    
    best = data(i+60,2);
    xneg = best - data(i+60,4);
    xpos = data(i+60,3) - best;
    
    plot(best,i-58,'or','markerfacecolor','r')
    
    
    errorbar(best,i-58,xneg,xpos,'horizontal','color','r')
    
end

axis([-0.1 5000000 0.3 2.2])


end