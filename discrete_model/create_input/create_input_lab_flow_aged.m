function create_input_lab_flow_aged
%% read data
% dir_data = 'procdata/05script/'

dir_data = '../../procdata/05script/';
dir_output = '../input/'

data = readtable([dir_data 'model_input_tompos_aged.csv']);


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
err_sd = zeros(n_clusters,n_utp);

data_norm = data_rid ./ repmat(str2double(data{:,end}),1,n_clusters) .* repmat(str2double(data{:,23}),1,n_clusters)



for j = 1:n_utp

    
    means(:,j) = mean(data_norm(time == unique_time(j),:),1)';
    sum(time == unique_time(j))
    err_sd(:,j) =std(data_norm(time == unique_time(j),:),1)/sqrt(sum(time == unique_time(j)));
   
end




dlmwrite([dir_output 'input_lab_num_tot_nodivide_aged.txt'],[unique_time,means',err_sd'])





end