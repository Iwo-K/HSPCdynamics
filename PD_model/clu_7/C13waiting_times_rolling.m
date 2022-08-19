function C13waiting_times_rolling

data_roll =  load('./PD_matrixes/results_rolling_clu_7.mat');
table_dpt = readtable('./tables/table_all_parameters_clu_7.csv');

dpt = sort(table_dpt{:,3});


AUC = 0;
solbest = 2;


vec = 0:1/298:1;

count = 0;

for int_start = dpt'
   
    int_stop = 1;
    
    
    while AUC<1
        
        count = count+1
        
        if count <=301
    yy = data_roll.sol{1,solbest}.x(count,:);
    
    
    
    AUC = trapz(vec(vec >= int_start & vec < int_stop),yy(vec >= int_start & vec < int_stop));
        else
            
            AUC = 1;
    
        end
    end
    
    dlmwrite('./tables/waiting_times_clu_7.txt',[int_start,count],'-append')

    
    count = count-1
    
    AUC = 0;
    
    
end



end

