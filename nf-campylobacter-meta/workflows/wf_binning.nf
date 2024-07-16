include { bwa } from  '../modules/mapping.nf'
include { metabat } from '../modules/binning.nf'

workflow METABAT_wf {
    take:
        ch_reads
        ch_ref
    main:
        bwa(ch_ref.map{it -> it[1]}, ch_reads)
        ch_metabat_input = bwa.out.reads.join(bwa.out.ref).join(bwa.out.bam)
        metabat(ch_metabat_input)
    emit:
        bin = metabat.out.bin
}