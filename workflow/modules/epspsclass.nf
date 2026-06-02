// workflow/modules/epspsclass.nf
// ================================
// Nextflow DSL2 module for EPSPS classification with EPSPSClass.
//
// Usage in your main workflow:
//
//   include { EPSPYCLASS_CLASSIFY } from './modules/epspsclass.nf'
//
//   workflow {
//       fasta_ch = Channel.fromPath("data/*.fasta")
//       EPSPYCLASS_CLASSIFY(fasta_ch)
//       EPSPYCLASS_CLASSIFY.out.tsv.view()
//   }

process EPSPYCLASS_CLASSIFY {
    tag "${meta.id}"
    label 'process_low'

    // Publish results to the output directory
    publishDir "${params.outdir}/classification", mode: 'copy'

    conda "conda-forge::biopython>=1.80"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_epspsclass.tsv"), emit: tsv
    path "${meta.id}_epspsclass.log",                   emit: log

    script:
    def threshold = params.epspsclass_threshold ?: 40.0
    """
    pip install epspsclass --quiet

    epspsclass classify \\
        --input  ${fasta} \\
        --output ${meta.id}_epspsclass.tsv \\
        --threshold ${threshold} \\
        --summary \\
        2> ${meta.id}_epspsclass.log
    """

    stub:
    """
    touch ${meta.id}_epspsclass.tsv
    touch ${meta.id}_epspsclass.log
    """
}


process EPSPYCLASS_CLASSIFY_STDIN {
    //
    // Streaming variant: reads from stdin, writes to stdout.
    // Use when chaining with other tools in a pipeline.
    //
    // Example:
    //   cat sequences.fasta | epspsclass classify -i - -o - | grep "^seq1"
    //
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/classification", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_epspsclass.tsv"), emit: tsv

    script:
    def threshold = params.epspsclass_threshold ?: 40.0
    """
    pip install epspsclass --quiet

    cat ${fasta} \\
        | epspsclass classify -i - -o - --threshold ${threshold} \\
        > ${meta.id}_epspsclass.tsv
    """
}
