#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -q fast.q
#$ -N sealion_test
#$ -pe smp 1

module load singularity/3.10.5_fix
CONTAINER="/share/scientific_bin/singularity/containers/SeaLion_container.sif"

cd code
singularity exec ${CONTAINER} sealion1
