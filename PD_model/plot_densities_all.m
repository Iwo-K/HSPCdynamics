function plot_densities_all

vec = 0:1/298:1;
n_grid = 300;
x = linspace(0,1,n_grid);
kk =0;
X = [vec';flipud(vec')]';

line0 = zeros(1,299);

start_i = [1 5 11 15 19 23 25 29 33];

stop_i = [4 10 14 18 22 24 28 32 36];

t = [3 7 12 27 49 76 112 161 269];

solbest=[ 5 2 2 2];

%%
figure(222)
clf

for time = 1:9
    kk=0;
    for clu = [ 6 10 7 11]       
        kk = kk+1;    
        
        data = load(['./clu_' num2str(clu) '/PD_matrixes/data_input_pd_clu_' num2str(clu) '.mat']);
        model =load(['./clu_' num2str(clu) '/PD_matrixes/results_clu_'  num2str(clu)  '.mat']);
        
        
        solbest_clu = solbest(kk);
        
        
        subplot(9,4, kk + 4*(time-1))

        
        u0 = [];
        jj = 0;
        for i = start_i(time):stop_i(time)
            jj = jj+1;
            temp = ksdensity(data.D.ind.hist{i},x,'support',[0,1],'function','pdf');
            u0(jj,:) = temp/trapz(x,temp);
        end
        
        
        u = mean(u0,1);
        u = 0.5*(u(1:n_grid-1)+u(2:n_grid));
        
        
        
        err = std(u0,1);
        err = 0.5*(err(1:n_grid-1)+err(2:n_grid));
        
   
        hold on
        
        
        Y = [u(1:299)'+err';flipud(u(1:299)'-err')]';
        
        fill(X,Y,'b','facecolor',[0 1 1]*0.8,'edgecolor',[0 0 1]*0.8)
        
        
        %         plot(vec,model.sol{1,solbest_clu}.p(time,:),'r','linewidth',2)
%         errorbar(vec,u,err,'b','linewidth',2)
%         plot(vec,u,'k','linewidth',2)
         plot(vec,model.sol{1,solbest_clu}.p(time,:),'r','linewidth',2)
        
        
        if kk == 1
            ylabel([num2str(t(time)) ' days'])
            
        end
        
        
        if time ==1
            title(['Cluster ' num2str(clu)])
        end
        
        
        if time  == 9
            xlabel('Pseudotime')
            
        end
        
    end
    
end

mkdir('figures')

set(gcf, 'PaperUnits', 'centimeters');
exportfig(gcf,'./figures/density_all.eps','FontMode', 'fixed','Fontsize',6,'color', 'cmyk','width',23,'height',30,'Renderer','painters','Lockaxes',0);%

end
