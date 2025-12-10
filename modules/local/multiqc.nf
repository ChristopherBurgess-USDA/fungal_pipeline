process MULTIQC {

    module 'multiqc'
    executor 'slurm'
    clusterOptions '--time=01:00:00 --cpus-per-task=2 --mem=4G -p ceres'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
        tuple val(seq_run), val(seq_direction)

    output:
        path("${seq_run}_${seq_direction}_multiqc_report.html")

    script:
    """
    multiqc ${params.outdir}/fastqc/${seq_direction}/ \
        --outdir . \
        --filename ${seq_run}_${seq_direction}_multiqc_report.html
    """
}