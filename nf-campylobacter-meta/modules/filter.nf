params.filter_reads = false

process filter_read {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 8
    
    conda 'bioconda::bbmap=38.90'

    input:
    tuple val(sample_id), path(reads), path(kraken2_output), val(species)

    output:
    tuple path("${sample_id}_${species}.tsv"), path("*.{fq.gz,fastq.gz,fasta}") optional true

    shell:
    if (params.filter_reads) {
        '''
        cat !{kraken2_output} | grep -i "!{species}" | cut -f 2 > campy
        if [[ -s campy ]]; then 
            repair.sh in=!{reads[0]} in2=!{reads[1]} out=R1.fq.gz out2=R2.fq.gz
            filterbyname.sh in=R1.fq.gz in2=R2.fq.gz names=campy out="!{sample_id}_!{species}_1.fq.gz" out2="!{sample_id}_!{species}_2.fq.gz" include=t qin=33
            mv campy !{sample_id}_!{species}.tsv
            rm R1.fq.gz R2.fq.gz
        fi
        '''
    } else {
        '''
        cat !{kraken2_output} | grep -i "!{species}" | cut -f 2 > campy
        if [[ -s campy ]]; then 
            filterbyname.sh in=!{reads} names=campy out="!{sample_id}_!{species}.fastq.gz" include=t qin=33
            mv campy !{sample_id}_!{species}.tsv
        fi
        '''
    }
    
}