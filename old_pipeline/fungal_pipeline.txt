
include { ITSXPRESS } from './modules/itsxpress.nf'
include { FASTQC } from './modules/fastqc.nf'
include { MULTIQC } from './modules/multiqc.nf'



workflow FUNGAL_PIPELINE {

    // Create a channel from the metadata CSV file with 4 columns
    metadata_ch = Channel
        .fromPath("${params.input_dir}/metadata.csv")
        .splitCsv(header:true)
        .map { row -> tuple(row.seq_run, row.sample_id, file("${params.input_dir}/${row.forward_reads}"), file("${params.input_dir}/${row.reverse_reads}")) }


    seq_direction_ch = Channel.of('R1','R2')

    seq_run_ch = Channel
        .fromPath("${params.input_dir}/metadata.csv")
        .splitCsv(header:true)
        .map { row -> row.seq_run }
        .distinct()
        .cross( seq_direction_ch )

    ITSXPRESS( metadata_ch )
    FASTQC( metadata_ch )
    MULTIQC( seq_run_ch )
}
