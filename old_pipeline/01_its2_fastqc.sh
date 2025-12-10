#!/bin/bash
#SBATCH --job-name="ky_its_qc"
#SBATCH --partition=ceres
#SBATCH --time=2:00:00
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 50 # number of logical cores/threads
#SBATCH --mail-user=christopher.burgess@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_out/%j-%x.out

module purge
module load fastqc multiqc
mkdir -p output_data/fastqc

for seq_run in run1 run2
do
    for fr in R1 R2
    do
        mkdir -p temp/01_fastqc/fastqc_${seq_run}/${fr}

        fastqc temp/00_itsxpress/reads/${seq_run}_*${fr}.fastq.gz \
            -f fastq \
            -t $SLURM_CPUS_ON_NODE \
            -o ./temp/01_fastqc/fastqc_${seq_run}/${fr}
        
        multiqc temp/01_fastqc/fastqc_${seq_run}/${fr}/ \
            --outdir temp/01_fastqc \
            --filename multiqc_${seq_run}_${fr}
    done
done

mv temp/01_fastqc/multiqc*.html output_data/fastqc