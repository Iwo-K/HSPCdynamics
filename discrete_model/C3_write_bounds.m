function C3_write_bounds

name_best = 'output/best.txt';


best=dlmread(name_best);

chibest = best(end);

for n=1:104
    
    
    fName=['./output/profile_likelihood/r_',num2str(n),'tot.txt'];
    
    
    if exist(fName,'file')==2
     
        
        data = dlmread(fName);
        
        
        data = sortrows(data,1);
        
        
        chi = data(:,end);
        
        data(chi > chibest + 3.8,:) ;
        
        data(chi > chibest + 3.8,:) = [];
        
        
        dlmwrite('./output/profile_likelihood/all_profiles.txt',data,'-append');
        
    end
    
    
    
    
end

delete './output/bounds.txt'

data = dlmread('./output/profile_likelihood/all_profiles.txt');




for lookAt = 1:104
    
    column = data(:,lookAt);
    
    
    
    [mini, index_mini] = min(column);
    [maxi, index_maxi] = max(column);
    
    
    
    chi_min = data(index_mini,end);
    chi_max = data(index_maxi,end);
    

    
    dlmwrite('./output/bounds.txt',[lookAt,best(lookAt),maxi,mini,chi_max,chi_min],'-append')
end

end
