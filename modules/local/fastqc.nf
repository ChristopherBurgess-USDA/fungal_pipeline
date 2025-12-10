process FASTQC {

    module 'fastqc'
    executor 'slurm'
    clusterOptions '--time=01:00:00 --cpus-per-task=2 --mem=4G -p ceres'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
        tuple val(seq_run), val(sample_id), path(forward_reads), path(reverse_reads)

    output:
        tuple val(sample_id), path("${seq_run}/${sample_id}_fastqc_R1.zip"), path("${seq_run}/${sample_id}_fastqc_R2.zip"), path("${seq_run}/${sample_id}_fastqc_R1.html"), path("${seq_run}/${sample_id}_fastqc_R2.html")

    script:
    """
    mkdir -p ${seq_run}/R1
    mkdir -p ${seq_run}/R2
    fastqc -f fastq --threads 2 --outdir . ${forward_reads} ${reverse_reads}

    mv ${forward_reads.baseName}_fastqc.zip ${seq_run}/R1/${sample_id}_fastqc_R1.zip
    mv ${reverse_reads.baseName}_fastqc.zip ${seq_run}/R2/${sample_id}_fastqc_R2.zip
    mv ${forward_reads.baseName}_fastqc.html ${seq_run}/R1/${sample_id}_fastqc_R1.html
    mv ${reverse_reads.baseName}_fastqc.html ${seq_run}/R2/${sample_id}_fastqc_R2.html
    """
}
