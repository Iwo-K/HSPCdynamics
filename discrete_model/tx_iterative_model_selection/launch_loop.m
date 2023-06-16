function launch_loop(launch)


vc = [ 0 1 2 3 4 5 6 7 8 9 10 11 12  16 24 26]; % cluster number that may change


%%

cluster_names = [0:12,14,16,20,24,25,26,28];

%% create differentiation matrix


M(1,:) = [0,4];
M(2,:) = [0,  5];
M(3,:) = [0,  8];
M(5,:) = [0, 22];


% from cluster 1

M(6,:) = [  1,  9];
M(7,:) = [  1, 12];


% from cluster 2

M(8,:) = [  2, 3];
M(9,:) = [  2,4];
M(10,:) = [  2, 5];
M(11,:) = [  2, 6];
M(12,:) = [  2, 12];
M(13,:) = [  2, 16];


% from cluster 3

M(14,:) = [  3, 10];
M(15,:) = [  3, 25];


% from cluster 4

M(16,:) = [4,  2];
M(17,:) = [4,  5];
M(18,:) = [4,  8];
M(19,:) = [4,  12];


% from cluster 5

M(20,:) = [  5, 2];
M(21,:) = [  5,4];
M(22,:) = [  5, 16];



% from cluster 6

M(23,:) = [  6, 24];


% from cluster 8

M(24,:) = [  8, 1];
M(25,:) = [  8,4];
M(26,:) = [  8, 7];
M(27,:) = [  8, 12];

% from cluster 9

M(28,:) = [  9, 11];

% from cluster 10

M(29,:) = [  10, 25];

% from cluster 11

M(30,:) = [  11, 20];

% from cluster 12

M(31,:) = [  12, 25];
M(32,:) = [  12, 26];

% from cluster 24

M(33,:) = [  24, 14];

% from cluster 16

M(34,:) = [  16, 14];
M(35,:) = [  16, 28];



% from cluster 30

M(36,:) = [21,0];
M(37,:) = [21,22];


% from cluster 40

M(4,:) = [22,  8];


%%

for j = 1+1000*(launch-1):1000*launch
    change_pop = [];

    pop = dec2bin(j)

    for mm = length(pop):-1:1

        if str2num( pop(mm))

            change_pop = [change_pop,vc(length(vc)+mm-length(pop))];

        end


    end




    para_change = [];

    for kk = change_pop

        para_change = [para_change;find(M(:,2) == kk)];
        para_change = [para_change;find(M(:,1) == kk)];

        para_change = [para_change;find(cluster_names==kk)+37];

        if kk == 0

            para_change = [para_change;58;4;37;59;60];

        end

    end
    para_change








    vec_para = para_change;
    vec_para = unique(vec_para);

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

    name_best = ['./bests/best_' num_hexa_dif '_' num_hexa_prol '_' num_hexa_log '.txt'];


    %%
    main_models(vec_para(:)',0);



    data2 = dlmread(name_best);
    chibest = data2(end);


    dlmwrite(['track.txt'],[change_pop,length(vec_para),chibest],'-append');


end
end


