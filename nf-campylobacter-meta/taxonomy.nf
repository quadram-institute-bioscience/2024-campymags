#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

//import modules

include {kraken2} from './modules/classifiers.nf'
include {krona} from './modules/classifiers.nf'
include {kraken2 as kraken2_standard} from './modules/classifiers.nf'
include {krona as krona_standard} from './modules/classifiers.nf'
include {bracken} from './modules/classifiers.nf' addParams(bracken_read_length: 150, taxlevel: 'S')
include {MetaPhlAn3} from './modules/classifiers.nf'

include {fastp;multiqc;checkm;quast;sendsketch} from './modules/qc.nf'
include {checkm as checkm_all} from './modules/qc.nf'


include {FILTER_SPECIES} from './workflows/wf_filter_campy.nf'
include {RELATIVE_ABUNDANCE} from './workflows/wf_abundance.nf'
include {REMOVED_HOST} from './workflows/wf_remove_host.nf'

if (params.single_end) {
    if (params.concatenate) {
        ch_reads = Channel.fromPath(params.reads, type: 'dir', checkIfExists: true)
                          .map {it -> tuple(it.getSimpleName(), it)}
    } else {
        ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
                          .map {it -> tuple(it.getSimpleName(), file(it))}
    }
} else {
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
}


kraken2_db = Channel.fromPath(params.kraken2_db)

ch_tax_level = Channel.of(params.tax_level.tokenize(',')).flatten()
bracken_db = Channel.fromPath(params.bracken_db)
ch_bracken = ch_tax_level.combine(bracken_db)


log.info 'TAXONOMY'
log.info "Input: ${params.contigs}"
log.info "single_end: ${params.single_end}"
log.info "DB location: ${params.kraken2_db}"
log.info "Out dir: ${params.outdir}"


process CONCATENATE {
    
    tag {dir_name}
    
    cpus 4

    input:
        tuple val(dir_name), path(dir)
    output:
        tuple val(dir_name), path("${dir_name}.fastq.gz")
    script:
    """
    zcat ${dir}/*.{fastq.gz,fq.gz,.fastq} | pigz - > ${dir_name}.fastq.gz
    """
}

process FASTQ_MERGE_LANES {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: 'copy'

    tag { sampleName }

    input:
    tuple val(sampleName), path(forward), path(reverse)
    
    output:
    tuple val(sampleName), path("${sampleName}_R1.fastq.gz"), path("${sampleName}_R2.fastq.gz")
    
    script:
    """
    zcat $forward | pigz -p 8 > ${sampleName}_R1.fastq.gz
    zcat $reverse | pigz -p 8 > ${sampleName}_R2.fastq.gz
    """
}

workflow MERGE_FASTQ_FOUR_LANES {
    take:
        ch_filePairs
    
    main:
        // Group 4 lane elements into 1 element in the channel
        ch_filePairs
                .map {
                it -> [it[0].replaceAll(~/\_L[1,2,3,4]/,"").replaceAll(~/${params.pattern}/,""), it[1][0], it[1][1]]
                }
                .groupTuple(by:0)
                .set { ch_reads_four_lanes }
        FASTQ_MERGE_LANES(ch_reads_four_lanes)
        ch_out = FASTQ_MERGE_LANES.out.map { it -> tuple(it[0], [it[1], it[2]]) }
    emit: 
        fastq = ch_out
}


workflow WIMP {
    take:
        ch_reads
    main:
        kraken2(ch_reads.combine(kraken2_db))
        ch_kraken2_output = kraken2.out.kraken2_output.map{it -> it[1]}
        if (params.krona_all) {
            ch_krona = ch_kraken2_output.collect()
        } else {
            ch_krona = ch_kraken2_output
        }
            krona(ch_krona)
}

workflow {
    main:
    if (params.concatenate && params.merge_lanes) { 
        log.info("Only accept concatenate or merge_lanes, not both")
        } 
    else {
        if (params.concatenate) {
            ch_reads1 = CONCATENATE(ch_reads)
        } else if (params.merge_lanes) {
            ch_reads1 = MERGE_FASTQ_FOUR_LANES(ch_reads)
        } else {
            ch_reads1 = ch_reads
        }
        if (params.abundance) {
            RELATIVE_ABUNDANCE(ch_reads1, kraken2_db, ch_bracken)
        }
        if (params.wimp) {
            WIMP(ch_reads1)
        }
        if (params.filter) {
            ch_kraken2_gtdb = Channel.fromPath(params.kraken2_db)
            ch_species = Channel.of(params.species.tokenize(',')).flatten()

            FILTER_SPECIES(ch_reads1, ch_species, ch_kraken2_gtdb)
        }
        
        if (params.remove_host) {
            ch_kraken2_host_db = Channel.fromPath("/dev/shm/kraken2_human_db_no_mask")
            REMOVED_HOST(ch_reads1, ch_kraken2_host_db)
        }
        if (params.metaphlan) {
            MetaPhlAn3(ch_reads1)
        }
    }
}
