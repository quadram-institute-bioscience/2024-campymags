include {kraken2} from '../modules/classifiers.nf' addParams(kraken2_use_name: true)
include {filter_read} from '../modules/filter.nf'

workflow FILTER_SPECIES {
    take: 
        ch_reads
        ch_species
        kraken2_db
    main:
        kraken2(ch_reads.combine(kraken2_db))
        filter_read(kraken2.out.reads.join(kraken2.out.kraken2_raw_output).combine(ch_species))
}