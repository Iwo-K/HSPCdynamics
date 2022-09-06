#!/bin/bash
#SBATCH -p icelake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --time 12:00:00
#SBATCH --job-name test3
#SBATCH --output test3_%J.log

####################################################################
# Modify this variables by adding the path to your singularity image
container="PATH TO THE CONTAINER"
container="/home/idk25/ik_rds/containers/rpy_v4/rpy_v4_p3_fix2.sif"

####################################################################

echo "Started:"; date
singularity exec ${container} R -e "knitr::spin('01_SS2_processing.R')"

pynotebooks="02_10x_processing
 04_integration
 05_label_discrete
 06_wot_cellrank
 07_projections
 22_Figure2_main
 23_Figure2_projections
 24_Figure3_main
 26_FigurePD
 36_3dhist
"

for i in $pynotebooks; do
    echo "Running: " $i
    singularity exec ${container} /usr/local/bin/jupytext --to notebook ${i}".py"
    singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace  \
                                    --execute --ExecutePreprocessor.timeout=0 ${i}".ipynb"
done

Rnotebooks="09_DE"
for i in $Rnotebooks
do
    echo "Running" $i
    singularity exec ${container} /usr/local/bin/jupytext --to notebook ${i}".R"
    singularity exec ${container} /usr/local/bin/jupyter nbconvert --to notebook --inplace  \
                                    --execute --ExecutePreprocessor.timeout=0 ${i}".ipynb"
done

echo "Finished:"; date
