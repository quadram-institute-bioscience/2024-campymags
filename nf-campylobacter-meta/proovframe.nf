#!/usr/bin/env nextflow

params.input = null // Should be  a CSV file.


process PROOVFRAME_MAP {
    
    tag {sample_id}
    
    conda 'bioconda::proovframe=0.9.7'
    
    cpus 64

    input:
    tuple val(sample_id), path(fasta), path(protein_db)
    output:
    tuple val(sample_id), path("${sample_id}_map.tsv")

    script:
    """
    proovframe-map -t ${task.cpus} -d ${protein_db} -o ${sample_id}_map.tsv ${fasta}
    """
}

process PROOVFRAME_FIX {
    publishDir "${params.outdir}", mode: 'copy'
    
    conda 'bioconda::proovframe=0.9.7'
    
    tag {sample_id}

    cpus 2

    input:
    tuple val(sample_id), path(fasta), path("${sample_id}_map.tsv")
    output:
    tuple val(sample_id), path("${sample_id}.frameshift.fixed.fasta")

    script:
    """
    proovframe-fix -o ${sample_id}.frameshift.fixed.fasta ${fasta} ${sample_id}_map.tsv
    """
}

ch_fasta = Channel.fromPath(params.input, checkIfExists: true)
                  .splitCsv(header: true)
                  .map {row -> {
                        [row.sample_id, row.fasta]
                    } 
                  }

ch_database = Channel.fromPath(params.database, checkIfExists: true)

ch_input = ch_fasta.combine(ch_database)

workflow {
    PROOVFRAME_MAP(ch_input)
    ch_in_proovframe_fix = ch_fasta.join(PROOVFRAME_MAP.out)
    PROOVFRAME_FIX(ch_in_proovframe_fix)
}