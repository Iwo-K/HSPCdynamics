function M = create_differentiation_matrix(cluster_names,n_clu)


% from cluster 0

M(1,:) = [find(cluster_names == 0),find(cluster_names == 4)];
M(2,:) = [find(cluster_names == 0),find(cluster_names == 5)];
M(3,:) = [find(cluster_names == 0),find(cluster_names == 8)];
M(5,:) = [find(cluster_names == 0), n_clu + 2];


% from cluster 1

M(6,:) = [find(cluster_names == 1),find(cluster_names == 9)];
M(7,:) = [find(cluster_names == 1),find(cluster_names == 12)];


% from cluster 2

M(8,:) = [find(cluster_names == 2),find(cluster_names == 3)];
M(9,:) = [find(cluster_names == 2),find(cluster_names == 4)];
M(10,:) = [find(cluster_names == 2),find(cluster_names == 5)];
M(11,:) = [find(cluster_names == 2),find(cluster_names == 6)];
M(12,:) = [find(cluster_names == 2),find(cluster_names == 12)];
M(13,:) = [find(cluster_names == 2),find(cluster_names == 16)];


% from cluster 3

M(14,:) = [find(cluster_names == 3),find(cluster_names == 10)];
M(15,:) = [find(cluster_names == 3),find(cluster_names == 25)];


% from cluster 4

M(16,:) = [find(cluster_names == 4),find(cluster_names == 2)];
M(17,:) = [find(cluster_names == 4),find(cluster_names == 5)];
M(18,:) = [find(cluster_names == 4),find(cluster_names == 8)];
M(19,:) = [find(cluster_names == 4),find(cluster_names == 12)];


% from cluster 5

M(20,:) = [find(cluster_names == 5),find(cluster_names == 2)];
M(21,:) = [find(cluster_names == 5),find(cluster_names == 4)];
M(22,:) = [find(cluster_names == 5),find(cluster_names == 16)];



% from cluster 6

M(23,:) = [find(cluster_names == 6),find(cluster_names == 24)];


% from cluster 8

M(24,:) = [find(cluster_names == 8),find(cluster_names == 1)];
M(25,:) = [find(cluster_names == 8),find(cluster_names == 4)];
M(26,:) = [find(cluster_names == 8),find(cluster_names == 7)];
M(27,:) = [find(cluster_names == 8),find(cluster_names == 12)];

% from cluster 9

M(28,:) = [find(cluster_names == 9),find(cluster_names == 11)];

% from cluster 10

M(29,:) = [find(cluster_names == 10),find(cluster_names == 25)];

% from cluster 11

M(30,:) = [find(cluster_names == 11),find(cluster_names == 20)];

% from cluster 12

M(31,:) = [find(cluster_names == 12),find(cluster_names == 25)];
M(32,:) = [find(cluster_names == 12),find(cluster_names == 26)];



% from cluster 16

M(34,:) = [find(cluster_names == 16),find(cluster_names == 14)];
M(35,:) = [find(cluster_names == 16),find(cluster_names == 28)];

% from cluster 24

M(33,:) = [find(cluster_names == 24),find(cluster_names == 14)];

% from cluster 30

M(36,:) = [n_clu + 1,find(cluster_names == 0)];
M(37,:) = [n_clu + 1,n_clu + 2];


% from cluster 40

M(4,:) = [n_clu + 2,find(cluster_names == 8)];


end