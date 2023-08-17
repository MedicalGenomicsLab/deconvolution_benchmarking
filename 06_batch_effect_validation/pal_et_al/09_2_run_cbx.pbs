#!/bin/bash
#######################
#PBS -N CBX-9
#PBS -l walltime=72:00:00
#PBS -l ncpus=2,mem=50Gb 
#PBS -j oe
#PBS -m ae
#PBS -M Khoa.TranQIMRBerghofer.edu.au
#PBS -W umask=0027
#PBS -J 1-19
#######################

module load singularity/3.7.1

PURITY_LEVELS=(
        0.05
        0.1
        0.15
        0.2
        0.25
        0.3
        0.35
        0.4
        0.45
        0.5
        0.55
        0.6
        0.65
        0.7
        0.75
        0.8
        0.85
        0.9
        0.95
)
PUR_LVL="${PURITY_LEVELS[$PBS_ARRAY_INDEX-1]}"

# Set working directory
WORK_DIR="???/deconvolution_benchmarking/06_batch_effect_validation/pal_et_al/data/cbx"
cd $WORK_DIR

# Specify paths to data files and output
SC_REF_PATH="scRNA_ref.tsv"
GENE_COUNTS_FILE="test_counts_${PUR_LVL}_pur_lvl.tsv"
OUT_DIR="./results/${PUR_LVL}"
mkdir -p $OUT_DIR

# Run CBX
singularity run \
-B ./:/src/data \
-B $OUT_DIR:/src/outdir \
--pwd /src \
???/fractions_latest.sif \
--username <EMAIL> \
--token <TOKEN> \
--single_cell TRUE \
--rmbatchSmode FALSE \
--perm 100 \
--refsample $SC_REF_PATH \
--mixture $GENE_COUNTS_FILE
