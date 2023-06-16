function analyse_track


data = dlmread('track.txt');


data = unique(data,"rows");

for j = 1: size(data,1)
    
    line = data(j,:);
    
    ii = find(line>0,1, 'last');
    
    if ii < length(line)
        
        line(ii+1:end) = [];
    end
    
    
    np(j) = length(line)-2;
    k(j) = line(end-1);
    chi(j) = line(end);
    
    a(j) = akaike(chi(j),63,k(j));
    
    
    
end

M = [np',k',chi',a'];

M(:,end+1) = (1:size(M,1))';




%%
% M(:,1:end-1)
% M = M( (M(:,4)<=115.8),:)
 M = M( (M(:,1)==5)  &  (M(:,3)<=70) &  (M(:,4)<=169) ,:)

 % M = M( (M(:,1)==5)  &  (M(:,3)<=70) ,:)
data(M(:,end),:)
min(M(:,3))






for i = 1:2%:size(M,1)
    
    line = data(M(i,end),:);
    
    ii = find(line>0,1, 'last');
    
    if ii < length(line)
        
        line(ii+1:end) = [];
    end
    
    
    vec_pop = line(1:end-2);
    
    vec_para = convert_pop_para(vec_pop);
    
return

    main_models(vec_para,1)
    
   
    
end



end