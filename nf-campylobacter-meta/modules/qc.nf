process fastp {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.{gz,json}",mode: "copy" 
    
    label 'high'

    tag {sample_id}
    
    conda 'bioconda::fastp=0.22.0'

    cpus 16
    memory '64.GB'

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastp.fq.gz"), path("${sample_id}_R2.fastp.fq.gz"), emit: reads
    path "${sample_id}.fastp.json", emit: json
    
    script:
    """
    fastp -5 -3 -r -M 20 --cut_front_window_size=1 --cut_tail_window_size=1 --cut_right_window_size=5 -w ${task.cpus} -z 6 \
    -i ${reads[0]} -I ${reads[1]} \
    -o ${sample_id}_R1.fastp.fq.gz -O ${sample_id}_R2.fastp.fq.gz \
    -j ${sample_id}.fastp.json
    """
}

process porechop {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::porechop=0.2.4'

    cpus 16
    memory '64.GB'
    
    input:
    tuple val(sample_id), path(fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_porechop.fastq.gz"), emit: trimmed
    path("${sample_id}_porechop.log.gz"), emit: log
    script:
    kit_name = params.kit_name ? "${params.kit_name}" : "auto"

    """
    porechop \
    -i ${fastq} \
    --kit_name ${kit_name} \
    --discard_middle \
    -t ${task.cpus} \
    -o '${sample_id}_porechop.fastq.gz'
    > ${sample_id}_porechop.log
    gzip ${sample_id}_porechop.log
    """
}



    // --detect_adapter_for_pe \
process filtlong {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.gz",mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::filtlong=0.2.0'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtlong.fastq.gz")
    
    script:
    """
    filtlong --min_length $params.min_length --keep_percent 95 ${reads} | gzip > ${sample_id}_filtlong.fastq.gz 
    """
}


process nanoq {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.gz",mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::nanoq=0.1.0'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_fil.fastq.gz")
    
    script:
    //Assume inputs are gzip
    """
    zcat ${reads} | nanoq -l ${params.min_length} -q ${params.min_qscore} | pigz - > ${sample_id}_fil.fastq.gz
    """
}


process multiqc {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.html",mode: "copy"
    
    tag "Multiqc"
    
    conda 'bioconda::multiqc=1.15'

    errorStrategy 'ignore'

    input:
    path(input)

    output:
    path "multiqc_report.html" optional true

    script:
    """
    multiqc --interactive .  
    """
}

process checkm {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    conda 'bioconda::checkm-genome=1.0.11'
    
    cpus 16
    memory '480 GB'
    
    input:
    path(bins)
    
    output:
    path("checkm_output"), emit: checkm
    path("*.tsv"), emit: report

    script:
    """
    checkm lineage_wf -t ${task.cpus} -x ${params.contigs_extension} --tab_table . checkm_output > checkm_output.tsv
    grep -v -e "\\[" checkm_output.tsv | csvtk -t cut -f 1,2,3,12-14 > checkm_output_mqc.tsv
    """
}

process busco {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: {filename -> if (filename.endsWith(".tsv")) { 
        "./reports/${filename}"
    } else { 
        "./${filename}"
        }
    }
    
    errorStrategy 'ignore'

    tag {sample_id}
    conda 'bioconda::busco=3.0.2'
    
    cpus 8
    
    memory '32.GB'
    
    input:
    tuple val(sample_id), path(contigs)
    
    output:
    path("${sample_id}")
    path("${sample_id}_busco.tsv"), emit: report
    path("${sample_id}/short_summary.*"), emit: summary
    
    script:
    offline = params.busco_offline ? "--offline" : ""
    lineage = params.busco_lineage ? "--lineage ${params.busco_lineage}" : "--auto-lineage-prok"

    """
    export TMPDIR=.;
    busco \
    -i ${contigs} \
    -o ${sample_id} \
     ${lineage} \
    -m geno \
    -c ${task.cpus} \
    ${offline} \
    --download_path ${params.busco_db_location}
    grep -e "C:" ${sample_id}/short_summary.specific* > ${sample_id}_busco
    awk '{print \$0="${sample_id}"\$0}' ${sample_id}_busco > ${sample_id}_busco.tsv
    """
}


process quast {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {"running"}
    
    conda 'bioconda::quast=5.0.2'

    cpus 8

    input:
    path(contigs)
    output:
    path("quast_output")
    path("quast_output/report.tsv"), emit: report
    
    script:
    """
    quast.py --min-contig 100 -t ${task.cpus} -o quast_output ${contigs}
    """
}

process sendsketch {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}

    cpus 4

    memory '16 GB'

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}.tsv")
    
    script:
    """
    export _JAVA_OPTIONS="-Xmx1g -Xms512m"
    sendsketch.sh in=${contigs} out=result
    tail -n +4 result | head -n 1 | cut -f 12 > ${sample_id}.tsv
    """
}

process gunc {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12
    memory '200.GB'

    input:
    tuple val(sample_id), path("contigs.fa")
    
    output:
    path("${sample_id}")
    script:
    """
    mkdir ${sample_id}
    ln -sf contigs.fa ${sample_id}.fa
    gunc run --threads ${task.cpus} \
    -i ${sample_id}.fa \
    --use_species_level \
    --out_dir ${sample_id} \
    --detailed_output \
    --sensitive \
    -r ${params.gunc_db}
    """
}


process rasusa {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag {sample_id}
    
    cpus 8
    memory '32.GB'

    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}_subsample.fastq.gz")
    script:
    """
    rasusa \
    --coverage ${params.coverage} \
    --genome-size ${params.genome_size} \
    --input ${fastq} \
    --output-type g \
    --seed 2021 \
    --output ${sample_id}_subsample.fastq.gz
    """
}

process DECONTAMINATE {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}",mode: "copy" 
    
    tag {sample_id}
    
    cpus 8

    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz")
    shell:
    '''
    repair.sh in1=!{reads[0]} in2=!{reads[1]} out1=R1.fq.gz out2=R2.fq.gz -da
    bowtie2 --threads !{task.cpus} --sensitive-local -x !{params.host_ref} -1 R1.fq.gz -2 R2.fq.gz -S temp.sam
    EXITCODE=$?
    if [ $EXITCODE -ne 0 ]; then
        ln -fs !{reads[0]} !{sample_id}_R1.fastq.gz
        ln -fs !{reads[1]} !{sample_id}_R2.fastq.gz
    else
        samtools view -@ !{task.cpus} -Sb temp.sam | samtools fastq -@ !{task.cpus} -f 4 -1 !{sample_id}_R1.fastq.gz -2 !{sample_id}_R2.fastq.gz
        rm temp.sam
    fi
    '''
}