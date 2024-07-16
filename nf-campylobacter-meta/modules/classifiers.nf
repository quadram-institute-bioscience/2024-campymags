process remove_host {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.fastq.gz",mode: "copy"
  label 'kraken2'
  
  tag {sample_id}
  
  conda 'bioconda::kraken2=2.1.1'

  cpus 8

  input:
  tuple val(sample_id), path(reads), path(db)

  output:
  tuple val(sample_id), path("*.fq.gz"), emit: host_fastq
  tuple val(sample_id), path("*.fastq.gz"), emit: non_host_fastq
  path "${sample_id}.kraken2.report", emit: report
  script:
  pe = params.single_end ? "" : "--paired"
  classified = params.single_end ? "${sample_id}.host.fq" : "${sample_id}.host#.fq"
  unclassified = params.single_end ? "${sample_id}.fastq" : "${sample_id}#.fastq"
  
  """
    hostname > hostname
    kraken2 \\
        --db $db \\
        --confidence $params.confidence \\
        --threads $task.cpus \\
        --unclassified-out $unclassified \\
        --classified-out $classified \\
        --report ${sample_id}.kraken2.report \\
        --report-zero-counts \\
        $pe \\
        --gzip-compressed \\
        $reads > /dev/null
    pigz -p $task.cpus *.fastq
    pigz -p $task.cpus *.fq
  """  
}


process kraken2 {
    publishDir "${params.outdir}", mode: "copy",
     saveAs: {filename -> if ((filename =~ /_[1,2]\.fastq\.gz/) || (filename =~ /host.fastq\.gz/) || filename.endsWith(".report") || filename.endsWith(".output"))
                            { if (params.remove_host) {
                              "removed_host/${filename}"
                              } else {
                                if (filename.endsWith(".report")) {
                                  "${task.process.replaceAll(":","_")}/${filename}"
                                } else if (filename.endsWith(".kraken2")) {
                                  "${task.process.replaceAll(":","_")}/${sample_id}.kraken2"
                                } else if (filename.endsWith(".output")) {
                                  "${task.process.replaceAll(":","_")}/${sample_id}.uncut.kraken2"
                                } 
                              }
                            }
                          else null}
    
    errorStrategy 'ignore'
    // maxForks 1

    label 'kraken2'

    tag {"${action}-->${sample_id}"}
    
    conda 'bioconda::kraken2=2.1.1'
    
    cpus 8
    
    memory '200 GB'

    input:
    tuple val(sample_id), path(reads), path(db)
    output:
    tuple val(sample_id), path(reads), emit: reads
    tuple val(sample_id), path("*.fq.gz"), emit: host_fastq  optional true
    tuple val(sample_id), path("*.fastq.gz"), emit: non_host_fastq optional true
    tuple val(sample_id), path("${sample_id}.kraken2.report"), emit: kraken2_report
    tuple val(sample_id), path("${sample_id}.kraken2") , emit: kraken2_output optional true
    tuple val(sample_id), path("kraken2.output"), emit: kraken2_raw_output optional true

    script:
    action = params.remove_host ? "removing host" : "classifying"
    pe = params.single_end ? "" : "--paired"
    classified = params.single_end ? "${sample_id}.host.fq" : "${sample_id}.host#.fq"
    unclassified = params.single_end ? "${sample_id}.non_host.fastq" : "${sample_id}#.non_host.fastq"
    write_fastq = params.remove_host ? "--unclassified-out $unclassified --classified-out $classified" : ""
    compress_output = params.remove_host ? "pigz -p $task.cpus *.fastq; pigz -p $task.cpus *.fq" : ""
    kraken2_output = (params.remove_host || params.kraken_output) ? "/dev/null" : "kraken2.output"
    krona_output = params.remove_host ? "" : "cat kraken2.output | cut -f 2,3 > ${sample_id}.kraken2"
    use_name = params.kraken2_use_name ? "--use-names" : ""
    """
    kraken2 \\
        --db $db \\
        --confidence $params.confidence \\
        --threads $task.cpus \\
        --memory-mapping \\
        --minimum-hit-groups 4 \\
        $write_fastq \\
        --report ${sample_id}.kraken2.report \\
        $use_name \\
        $pe \\
        $reads > ${kraken2_output}
    ${compress_output}
    ${krona_output}
    """
}


process bracken {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {"${sample_id}-${tax_level}"}

    conda 'bioconda::bracken=2.6.0'
    // label 'kraken2'
    
    input:
    tuple val(sample_id), path("${sample_id}_kraken2.report"), val(tax_level), path(bracken_db)
    output:
    tuple val(sample_id), path("${sample_id}_bracken_${tax_level}.tsv")
    
    script:
    """
    bracken \
        -d ${bracken_db} \
        -r ${params.bracken_read_length} \
        -i ${sample_id}_kraken2.report \
        -l ${tax_level} \
        -o ${sample_id}_bracken_${tax_level}.tsv \
        -w ${sample_id}_bracken.kraken2_report
    """
}

process bracken_sort {
    // publishDir 
    
    tag {sample_id}
    
    cpus 1

    input:
    tuple val(sample_id), path("${sample_id}_bracken.tsv")
    output:
    path("${sample_id}_bracken_sorted.tsv")

    script:
    """
    csvtk -t mutate2 -T -n sample -e "'${sample_id}'" "${sample_id}_bracken.tsv" \
    | csvtk -t sort -k fraction_total_reads:r \
    | csvtk head -n 3 > ${sample_id}_bracken_sorted.tsv

    """
}

process bracken_collect {
    publishDir "${params.outdir}", mode: "copy" 
    
    tag {"running"}
    
    cpus 1

    input:
    path bracken

    output:
    path("bracken_all_samples.tsv")

    script:
    """
    csvtk concat -t ${bracken} > bracken_all_samples.tsv
    """
}

process krona {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
  tag {"running"}
  
  memory {params.krona_all ? '164 GB' : ' 128 GB'}
  
  input:
  path(kraken2)
  output:
  path("*.html")
  
  script:
  if (params.krona_all) {

      """
      ktImportTaxonomy -tax $params.krona_taxonomy ${kraken2}
      """
  } else {
      """
      ktImportTaxonomy -tax $params.krona_taxonomy -o ${kraken2.baseName}.krona.html ${kraken2}
      """
  }
  
}


process MetaPhlAn3 {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 8

    input:
        tuple val(sample_id), path(reads)
    output:
        path("${sample_id}_profile.txt"), emit: profile
        path("${sample_id}.bowtie2out.txt"), emit: bowtie2output
    
    script:

    """
        metaphlan ${reads[0]},${reads[1]} \
        --input_type fastq \
        --bowtie2db ${params.bowtie2db} \
        --nproc $task.cpus \
        --sample_id ${sample_id} \
        --bowtie2out ${sample_id}.bowtie2out.txt \
        > ${sample_id}_profile.txt
    """
}


process MMSEQS_EASY_TAXONOMY {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12
    memory '256.G'

    input:
        tuple val(sample_id), path(fasta), path(mmseqs_db)
    output:
        tuple val(sample_id), path("${sample_id}_lca.tsv")               , emit: lca
        tuple val(sample_id), path("${sample_id}_report")                , emit: report
        tuple val(sample_id), path("${sample_id}_tophit_aln")            , emit: tophit_aln
        tuple val(sample_id), path("${sample_id}_tophit_report")         , emit: tophit_report
    script:
    def mmseqs_db_name = mmseqs_db.getBaseName()

    """
    mmseqs easy-taxonomy \
    ${fasta} ${mmseqs_db}/${mmseqs_db_name} ${sample_id} tmp --filter-hits 1 --threads ${task.cpus} --orf-filter 1 --sort-results 1 --lca-mode 4 --lca-ranks species,genus --tax-lineage 1
    """
}