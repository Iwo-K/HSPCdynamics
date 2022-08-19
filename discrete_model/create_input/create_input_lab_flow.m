function create_input_lab_flow_rel_to_day3
%% read data
% dir_data = 'procdata/05script/'

dir_data = '../../procdata/05script/';
dir_output = '../input/'


data = readtable([dir_data 'model_input_tompos.csv']);


remove = isnan(str2double(data{:,23}));
data(remove,:) = [];

data_rid = data{:,2:21};

n_clusters = size(data_rid,2);


time = str2double(data{:,25});



unique_time = unique(time);
n_tp = length(time);


n_utp = length(unique_time);



%% averages over time


means = zeros(n_clusters,n_utp);
err_pv = zeros(n_clusters,n_utp);
sum_iter = zeros(1,n_clusters);



data_norm = data_rid ./ repmat(str2double(data{:,end}),1,n_clusters) .* repmat(str2double(data{:,23}),1,n_clusters);



for j = 1:n_utp
   
sum(time == unique_time(j))
    means(:,j) = nanmean(data_norm(time == unique_time(j),:),1)';
    
    sum_iter = sum_iter + nanstd(data_norm(time == unique_time(j),:),1).^2*(sum(time == unique_time(j))-1);
    
end




pv_t = sum_iter/(n_tp-length(unique_time));



for j = 1:length(unique_time)
 
    
    err_pv(:,j) = sqrt(pv_t)*1./sqrt(sum(time == unique_time(j))-1);
end

dlmwrite([dir_output 'input_lab_num_tot_nodivide.txt'],[unique_time,means',err_pv'])





end