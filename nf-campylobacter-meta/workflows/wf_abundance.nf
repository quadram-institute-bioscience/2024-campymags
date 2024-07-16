
include {kraken2} from '../modules/classifiers.nf'
include {bracken} from '../modules/classifiers.nf'

workflow RELATIVE_ABUNDANCE {
    take:
        ch_reads
        ch_kraken2_db
        ch_bracken
    main:
        kraken2(ch_reads.combine(ch_kraken2_db))
        ch_input_bracken = kraken2.out.kraken2_report.combine(ch_bracken)
        bracken(ch_input_bracken)
}

