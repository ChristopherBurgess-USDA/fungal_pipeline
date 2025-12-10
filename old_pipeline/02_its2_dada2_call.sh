#!/bin/bash
#SBATCH --job-name="ky_its2_dada"
#SBATCH --partition=ceres
#SBATCH --time=2:00:00
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 72 # number of logical cores/threads
#SBATCH --mail-user=Bryan.Emmett@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_out/%j-%x.out

module purge
module load gcc glpk gmp r miniconda # For some reason need gcc glpk and gmp for phyloseq

mkdir -p temp/03_dada temp/05_funguild output_data/asv mlsh_import_data

#*****************************************************************************#
#                   Importing Data from MLSH sequencing                       #
#*****************************************************************************#

cp ../../MLSH/02_ITS2/temp/03_dada/04_seqtab*nochim_run2.rds mlsh_import_data
cp ../../MLSH/02_ITS2/temp/03_dada/04_seqtab*nochim_run3.rds mlsh_import_data
cp ../../MLSH/02_ITS2/output_data/ky21_mlsh_seq_run_meta_data.csv mlsh_import_data


Rscript --vanilla 02-1_its2_dada2.R

Rscript --vanilla 02-2_its2_phyloseq_gen.R
