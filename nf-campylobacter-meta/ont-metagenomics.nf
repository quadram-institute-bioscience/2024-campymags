#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

MODULES = './modules/'


include {flye as flye_meta} from MODULES + 'assembler.nf' addParams(flye_options: '--meta')
include {busco as busco_medaka} from MODULES + 'qc.nf'
include {medaka} from MODULES +  'polishing.nf' addParams(medaka_model: 'r941_prom_hac_g507') //r941_prom_hac_g507 r941_min_hac_g507
include {filtlong} from MODULES + 'qc.nf'  addParams(min_length: '800')
include {porechop} from MODULES + 'qc.nf'

ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
                    .map {it -> tuple(it.baseName.replaceAll(/.fastq/,''), file(it))}
                    .filter{!(it[0] =~ /unclassified/)}

ch_genome_size = Channel.of(params.genome_size)

kraken2_db = Channel.fromPath(params.kraken2_db)

log.info 'C A M P Y N A N O P O R E'
log.info "Input: ${params.reads}"
log.info "Out dir: ${params.outdir}"


workflow {

    filtlong(ch_reads)

    porechop(filtlong.out)

    flye_meta(porechop.out.trimmed)

    medaka(porechop.out.trimmed.join(flye_meta.out.contigs))

    busco_medaka(medaka.out)

}