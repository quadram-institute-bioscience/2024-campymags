params.pilon_output = false
params.memory = '32G'
params.iteration = 20

process medaka {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    label 'medaka'
    
    tag {sample_id}
    
    cpus 12

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_medaka.fasta")
    
    script:
    """
    medaka_consensus -i ${reads} -d ${contigs} -o output -t ${task.cpus} -m ${params.medaka_model}
    mv output/consensus.fasta ${sample_id}_medaka.fasta
    """
}

process pilon_raw {
    publishDir "${params.outdir}", mode: "copy",
    saveAs: {filename -> if (params.pilon_output) { 
                "${task.process.replaceAll(":","_")}/${filename}"
                } 
            else {
                null
                }
            }

    tag {sample_id}
    
    cpus 8
    // memory params.memory

    input:
    tuple val(sample_id), path(short_reads), path(contigs) 

    output:
    tuple val(sample_id), path("${sample_id}_pilon.fasta"), emit: contigs
    path("*.changes")

    script:
    """
    export _JAVA_OPTIONS="-Xmx${params.memory} -Xms512m"
    bwa index ${contigs}
    
    bwa mem -M -t ${task.cpus} ${contigs} ${short_reads} | samtools view -bS -| samtools sort > alignments.bam

    samtools index alignments.bam

    pilon --genome ${contigs} --frags alignments.bam --changes \
    --output ${sample_id}_pilon --fix all
    """
}


process pilon {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    path "*.fasta"
    path "*.changes", optional: true
    
    script:
    """
    pilon.py -t ${task.cpus} -n ${params.iteration} ${reads[0]} ${reads[1]} ${contigs}
    mv final.polished.fasta  ${sample_id}_pilon.fasta
    """
}

process POLYPOLISH {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag {sample_id}
    
    cpus 8
    memory '32.GB'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_polypolish.fasta")
    script:
    """
    bwa index ${contigs}
    bwa mem -t ${task.cpus} -a ${contigs} ${reads[0]} > alignments_1.sam
    bwa mem -t ${task.cpus} -a ${contigs} ${reads[1]} > alignments_2.sam
    polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
    polypolish ${contigs} filtered_1.sam filtered_2.sam > ${sample_id}_polypolish.fasta 
    """
}

process POLCA {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 12

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    path "*.PolcaCorrected.fa", emit: contigs
    path "*.report", emit: report
    script:
    """
    polca.sh -a ${contigs} -r \'${reads}\' -t ${task.cpus} -m 1G
    """
}