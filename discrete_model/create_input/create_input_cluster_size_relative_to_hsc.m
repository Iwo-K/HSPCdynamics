function create_input_cluster_size_relative_to_hsc
%% read data

dir_data = '../../procdata/05script/';
dir_output = '../input/'

data = readtable([dir_data 'model_input_tomneg.csv']);

clu_hsc = 0 + 1;
data_rid = data{2:end,2:21};



n_clusters = size(data_rid,2);




time = str2double(data{2:end,25});


unique_time = unique(time);
n_tp = length(time);


n_utp = length(unique_time);

%% averages over time


means = zeros(n_clusters,n_utp);
err_pv = zeros(n_clusters,n_utp);
sum_iter = zeros(1,n_clusters);


data_norm = data_rid./repmat(data_rid(:,clu_hsc),1,n_clusters);

for j = 1:n_utp
   
    
    
    means(:,j) = mean(data_norm(time == unique_time(j),:))';
    
    sum_iter = sum_iter + std(data_norm(time == unique_time(j),:)).^2*(sum(time == unique_time(j))-1);
    
end



pv_t = sum_iter/(n_tp-length(unique_time));



for j = 1:length(unique_time)
 
    
    err_pv(:,j) = sqrt(pv_t)*1./sqrt(sum(time == unique_time(j))-1);
end


dlmwrite([dir_output 'input_cluster_size_rel_hsc_over_time_pv.txt'],[unique_time,means',err_pv'])


end