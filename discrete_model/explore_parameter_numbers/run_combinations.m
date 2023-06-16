function run_combinations(num)

vecs = dlmread('./output/df22.csv');

size(vecs)

for i = num:100:12000
i

vec = vecs(i,2:end);

main_combinations(vec)

end


end
