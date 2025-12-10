
process ITSXPRESS {

    module 'itsxpress'
    executor 'slurm'
    clusterOptions '--time=02:00:00 --cpus-per-task=4 --mem=12G -p ceres'

    publishDir "${params.outdir}/itsxpress", mode: 'copy'

    input:
        tuple val(seq_run), val(sample_id), path(forward_reads), path(reverse_reads)

    output:
        tuple val(sample_id), path("${seq_run}_${sample_id}_itsxpress_R1.fastq.gz"), path("${seq_run}_${sample_id}_itsxpress_R2.fastq.gz"), path("reads_track.csv")

    script:
    """
    itsxpress \
        --fastq ${forward_reads} \
        --fastq2 ${reverse_reads} \
        --outfile ${seq_run}_${sample_id}_itsxpress_R1.fastq.gz \
        --outfile2 ${seq_run}_${sample_id}_itsxpress_R2.fastq.gz \
        --region ITS2 \
        --threads 4 \
        --log ${seq_run}_${sample_id}_itsxpress.log

    pipeline_id="${seq_run}_${sample_id}" # quick awk fix
    tail -5 ${seq_run}_${sample_id}_itsxpress.log |\
        head -4 |\
        awk '{printf ",%s", substr($NF, 1, length($NF)-1)}' |\
        awk -v x="$pipeline_id" '{print x $0}' >>\
        reads_track.csv
    """
}
