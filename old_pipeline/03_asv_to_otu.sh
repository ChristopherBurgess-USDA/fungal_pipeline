#!/bin/bash
#SBATCH --job-name="ky_its_vsearch"
#SBATCH -p mem-low # change queue based on availability 
#SBATCH --nodes=1   # number of nodes
#SBATCH --time=2:00:00
#SBATCH -n 60 # number of logical cores/threads
#SBATCH --mail-user=christopher.burgess@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_out/%j-%x.out

### File Prep
module load vsearch qiime2

mkdir -p temp/04_vsearch

for seq_type in forward merged
do
#*****************************************************************************#
#              Clustering dada2 asv to otus                                   #
#*****************************************************************************#
    vsearch --cluster_smallmem output_data/asv/its2_asv_${seq_type}_rep_set.fasta \
        --id 0.97 \
        --threads $SLURM_CPUS_ON_NODE \
        --strand plus \
        --qmask none \
        --centroids output_data/its2_otu_${seq_type}_rep_set.fasta \
        --relabel otu_ \
        --usersort

#*****************************************************************************#
#              Mapping asvs to otus                                           #
#*****************************************************************************#
    vsearch --usearch_global output_data/asv/its2_asv_${seq_type}_rep_set.fasta \
        --db output_data/its2_otu_${seq_type}_rep_set.fasta \
        --id 0.97 \
        --uc temp/04_vsearch/01_asv_otu_${seq_type}_map.uc \
        --strand plus \
        --threads $SLURM_CPUS_ON_NODE

#*****************************************************************************#
#              Using Qiime2 for taxanomy                                      #
#*****************************************************************************#
    qiime tools import \
        --input-path output_data/its2_otu_${seq_type}_rep_set.fasta \
        --output-path temp/04_vsearch/02_its2_otu_${seq_type}_rep_set.qza \
        --type 'FeatureData[Sequence]'

    qiime feature-classifier classify-sklearn \
        --i-reads temp/04_vsearch/02_its2_otu_${seq_type}_rep_set.qza \
        --i-classifier ../../databases/unite_ver10_dynamic_04.04.2024-Q2-2024.2.qza \
        --p-n-jobs $SLURM_CPUS_ON_NODE \
        --o-classification temp/04_vsearch/03_its2_otu_${seq_type}_taxanomy.qza

    qiime tools export \
        --input-path temp/04_vsearch/03_its2_otu_${seq_type}_taxanomy.qza \
        --output-path temp/04_vsearch

    mv temp/04_vsearch/taxonomy.tsv \
        temp/04_vsearch/05_its2_otu_${seq_type}_unite-04-04-2024_taxanomy.tsv
done

module purge
module load gcc glpk gmp r

Rscript --vanilla 03-1_otu_phyloseq_gen.R
