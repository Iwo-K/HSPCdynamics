function C3_write_bounds
vec_para = [1,6,7,8,9,10,11,12,13,16,17,18,19,20,21,23,24,25,27,31,32,33,39,40,42,50,54];

vec_para = unique(vec_para);


vec_para = sort(vec_para);



vec_dif = sort(vec_para(vec_para<=37));


vec_prol = sort(vec_para((vec_para>37) & (vec_para<=58))-37);


vec_log = sort(vec_para(vec_para>58)-58);


%%

num_dec_dif = sum(2.^(37-vec_dif));


num_hexa_dif = dec2hex(num_dec_dif);


num_dec_prol = sum(2.^(21-vec_prol));


num_hexa_prol = dec2hex(num_dec_prol);


num_dec_log = sum(2.^(2-vec_log));


num_hexa_log = dec2hex(num_dec_log);


%%


name_best = ['./bests/best_' num_hexa_dif '_' num_hexa_prol '_' num_hexa_log '.txt']

%%

best=dlmread(name_best);

chibest = best(end);

for n=1:27
    
    
    fName=['./pl_tm_1/r_',num2str(n),'tot.txt'];
    
    
    if exist(fName,'file')==2
     
        
        data = dlmread(fName);
        
        
        data = sortrows(data,1);
        
        
        chi = data(:,end);
        
        data(chi > chibest + 3.8,:) ;
        
        data(chi > chibest + 3.8,:) = [];
        
        
        dlmwrite('./pl_tm_1/all_profiles.txt',data,'-append');
        
    end
    
    
    
    
end

delete './pl_tm_1/bounds.txt'

data = dlmread('./pl_tm_1/all_profiles.txt');

best_base =  dlmread('./input_data/best_base.txt');
min(data(:,end))

for lookAt = 1:27
    change = 1;
    
    column = data(:,lookAt);
    
    
    
    [mini, index_mini] = min(column);
    [maxi, index_maxi] = max(column);
    
    ratio = mini/best_base(vec_para(lookAt)+44);
    
    chi_min = data(index_mini,end);
    chi_max = data(index_maxi,end);
    
    if (maxi >= 0.0001) & (mini <= 0.0001)
change = 0;
ratio = 0;
    end

if maxi <= -0.0001
change  = -1;
ratio = maxi/best_base(vec_para(lookAt)+44);
end

    
    dlmwrite('./pl_tm_1/bounds.txt',[vec_para(lookAt),best(lookAt),maxi,mini,chi_max,chi_min,change,ratio],'-append')
end

end