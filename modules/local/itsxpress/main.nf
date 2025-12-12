// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process ITSXPRESS {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    
    tuple val(meta), path(reads)

    output:
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    tuple val(meta), path("*.{fastq}"), emit: fastq
    // TODO nf-core: Update the command here to obtain the version number of the software used in this module
    // TODO nf-core: If multiple software packages are used in this module, all MUST be added here
    //               by copying the line below and replacing the current tool with the extra tool(s)
    tuple val("${task.process}"), val('itsxpress'), eval("itsxpress --version"), topic: versions, emit: versions_itsxpress

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // added in ITS2 region to modules.config
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_command = "--fastq" + reads.collect{ " ${it}" }.join(" --fastq2")
    def reads_out = reads.withIndex().collect { entry, index -> [ entry, "itsxpress_${prefix}_${index + 1}.${entry.extension}" ] }.collect{ _in_name, out_name -> out_name }
    def output_command = "--outfile" + reads_out.collect{ " ${it}" }.join(" --outfile2")
    
    """

    ## Need to check syntax for 
   
    itsxpress \\
        $args \\ 
        -threads $task.cpus \\
        $input_command \\
        $output_command \\
        --log ${prefix}_itsxpress.log
    """
}