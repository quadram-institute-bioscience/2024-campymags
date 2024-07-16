
include {kraken2 as remove_host} from '../modules/classifiers.nf' addParams(remove_host: true)

workflow REMOVED_HOST {
    take:
        ch_reads
        ch_kraken2_db
    main:
        remove_host(ch_reads.combine(ch_kraken2_db))
}
