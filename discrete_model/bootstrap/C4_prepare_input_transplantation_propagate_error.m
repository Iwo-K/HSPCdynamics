function prepare_input_transplantation_propagate_error
%% read data

input_dir = '../../procdata/07script/';

data_df = readtable([input_dir 'tx_cellnos_permouse.csv'],'ReadRowNames',true);
data  = data_df{3:end-1,1:end-3};
l_m = size(data,1);



n(1) = sum(contains(data_df.Properties.VariableNames,'1Day'));
n(2) = sum(contains(data_df.Properties.VariableNames,'3Day'));
n(3) = sum(contains(data_df.Properties.VariableNames,'5Day'));
n(4) = sum(contains(data_df.Properties.VariableNames,'7Day'));


leiden_clusters = data_df{1,1:end-3};

cluster_names = [0:12,14,16,20,24,25,26,28];
n_clusters = length(cluster_names);


%% prepare matrix

vec_all = repmat(cluster_names',4,1);

vec_time = [ones(n_clusters,1);ones(n_clusters,1)*3;ones(n_clusters,1)*5;ones(n_clusters,1)*7];

mat_all = [vec_time, vec_all];



vec_time = [ones(n(1),1);ones(n(2),1)*3;ones(n(3),1)*5;ones(n(4),1)*7];

mat = [vec_time,leiden_clusters'];




v = ismember(mat_all,mat,'rows');

data_all = zeros(l_m,n_clusters*4);

data_all(:,v) = data;



for tp = 1:4
    
index_day(:,tp) = sum(data_all(:,1+(tp-1)*20:20*tp),2);

end
index_day_bool = index_day>0;

n_mice_time = sum(index_day_bool);



%% propagate error

for tp = 1:4
mean_clusters(tp,:) = mean(data_all(index_day_bool(:,tp),1+(tp-1)*n_clusters:n_clusters*tp));


end


err_mean_clusters = sqrt(mean_clusters); % poissonian error


mean_hsc = repmat(mean_clusters(:,1),1,n_clusters);

err_mean_hsc = repmat(err_mean_clusters(:,1),1,n_clusters);



mean_normalized = mean_clusters./mean_hsc;

prop_err = mean_normalized .* sqrt((err_mean_clusters./mean_clusters).^2 + (err_mean_hsc./mean_hsc).^2);

prop_err = prop_err./repmat(sqrt(n_mice_time'-1),1,n_clusters);
%% write output

dlmwrite('input_tx_boot_err_prop.txt',prop_err)


end