function count_combinations

cluster_names = [0:12,14,16,20,24,25,26,28];

n_clu = length(cluster_names);

mat = create_differentiation_matrix(cluster_names,n_clu);

count = 0;


for i1 = [1,0]
    for i2 = [1,0]
        for i3 = [1,0]
            for i4 = [1,0]
                for i5 = [1,0]
                    for i6 = [1,0]
                        for i7 = [1,0]
                            for i8 = [1,0]
                                for i9 = [1,0]
                                    for i10 = [1,0]
                                        for i11 = [1,0]
                                            for i12 = [1,0]
                                                for i13 = [1,0]
                                                    for i14 = [1,0]
                                                        for i15 = [1,0]
                                                            for i16 = [1,0]
                                                                for i17 = [1,0]
                                                                    for i18 = [1,0]
                                                                        for i19 = [1,0]
                                                                            for i20 = [1,0]
                                                                                for i21 = [1,0]
                                                                                    for i22 = [1,0]
                                                                                        for i23 = [1,0]
                                                                                            for i24 = [1,0]
                                                                                                for i25 = [1,0]
                                                                                                    for i26 = [1,0]
                                                                                                        for i27 = [1,0]
                                                                                                            for i28 = [1,0]
                                                                                                                for i29 = [1,0]
                                                                                                                    for i30 = [1,0]
                                                                                                                        for i31 = [1,0]
                                                                                                                            for i32 = [1,0]
                                                                                                                                for i33 = [1,0]
                                                                                                                                    for i34 = [1,0]
                                                                                                                                        for i35 = [1,0]
                                                                                                                                            for i36 = 1
                                                                                                                                                for i37 = [1,0]




                                                                                                                                                    vec = [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,i28,i29,i30,i31,i32,i33,i34,i35,i36,i37];

                                                                                                                                                    %  vec = ([0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,1,0,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1]);

                                                                                                                                                    sum(vec)


                                                                                                                                                    mat_sub = mat(vec == 1,:);

                                                                                                                                                    % uv = unique(mat_sub(:,2));

                                                                                                                                                    ll = [];
                                                                                                                                                    llp = [1;22];

                                                                                                                                                    % uu = 21;

                                                                                                                                                    while length(llp) > length(ll)

                                                                                                                                                        ll = llp;

                                                                                                                                                        new = mat_sub(find(ismember(mat_sub(:,1),llp)),:);

                                                                                                                                                        %uv = unique(new(:,2));

                                                                                                                                                        llp = unique([llp;new(:,2)]);


                                                                                                                                                    end




                                                                                                                                                    if length(llp) == 21


                                                                                                                                                        count = count +1

                                                                                                                                                        dlmwrite('./output/model_vectors.txt',vec,'-append')
                                                                                                                                                    end


                                                                                                                                          
                                                                                                                                                end


                                                                                                                                            end

                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                end

                                                                                                                            end
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                        end
                                                                                                    end

                                                                                                end
                                                                                            end
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end

                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


end