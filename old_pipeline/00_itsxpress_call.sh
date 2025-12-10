#!/bin/bash
#SBATCH --job-name="ky_itsxpress"
#SBATCH -p mem768 # change queue based on availability 
#SBATCH --nodes=1   # number of nodes
#SBATCH -n 60 # number of logical cores/threads
#SBATCH --mem-per-cpu 12GB
#SBATCH --mail-user=christopher.burgess@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_out/%j-%x.out

module purge

module load miniconda
source activate /project/mmicrobe/itsxpress-env

mkdir -p temp/00_itsxpress/reads temp/00_itsxpress/logs

echo sample_id,raw_forward,raw_reverse,itsxpress_forward,itsxpress_reverse > temp/00_itsxpress/reads_track.csv

for seq_run in run1 run2
do

    while read -r sample_id raw_forward raw_reverse
    do
    pipeline_id="${seq_run}_${sample_id}"
        itsxpress \
            --fastq  ${raw_forward} \
            --fastq2 ${raw_reverse} \
            --region ITS2 \
            --outfile temp/00_itsxpress/reads/${pipeline_id}_R1.fastq.gz \
            --outfile2 temp/00_itsxpress/reads/${pipeline_id}_R2.fastq.gz \
            --threads $SLURM_CPUS_ON_NODE \
            --log temp/00_itsxpress/logs/${pipeline_id}.txt
        
        tail -5 temp/00_itsxpress/logs/${pipeline_id}.txt |\
            head -4 |\
            awk '{printf ",%s", substr($NF, 1, length($NF)-1)}' |\
            awk -v x="$pipeline_id" '{print x $0}' >>\
            temp/00_itsxpress/reads_track.csv

    done < ../00_Raw_Data/ITS2_manifest_${seq_run}.txt
done

awk -v FS=, '{print $1}' temp/00_itsxpress/reads_track.csv >\
    temp/00_itsxpress/sample_names.txt