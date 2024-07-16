params.prokka_options = ""

process prokka {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"

    tag {sample_id}
    
    conda 'bioconda::prokka=1.14.5'

    errorStrategy 'ignore'
    
    cpus 8

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}")
    path "${sample_id}/${sample_id}.gff" , emit: roary
    
    script:
    """
    prokka --kingdom Bacteria ${params.prokka_options} --outdir ${sample_id} --prefix ${sample_id} ${contigs}
    """
}


process roary {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::roary=3.13.0'

    errorStrategy 'ignore'
    
    cpus 16

    input:
    path(gffs)
    output:
    path("roary_results")

    script:
    """
    roary -p ${task.cpus} -e --mafft -f roary_results ${gffs}
    """
}

process mlst {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}

    conda 'bioconda::mlst=2.19.0'

    errorStrategy 'ignore'
    
    cpus 1

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}_mlst.tsv")

    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.tsv
    """
}


process staramr {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::staramr=0.8.0'

    errorStrategy 'ignore'
    
    cpus 1

    memory '64.GB'

    input:
    path(contigs)
    output:
    path("staramr")

    script:
    """
    export BLASTDB_LMDB_MAP_SIZE=1000000
    staramr search --output-dir staramr *.fa*
    """
}



process amrfinderplus {
    publishDir "${params.outdir}/amrfinderplus", mode: "copy"
    
    tag {sample_id}
    
    conda 'bioconda::ncbi-amrfinderplus=3.10.1'

    errorStrategy 'ignore'
    
    cpus 2

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("*.fna")           , emit: fasta
    path("${sample_id}.tsv"), emit: tsv

    script:
    """
    amrfinder -t ${task.cpus} \
    --organism ${params.amrfinderplus_species} \
    --nucleotide ${contigs} \
    --report_common \
    -o ${sample_id}.tsv \
    --name ${sample_id} \
    --nucleotide_output ${sample_id}.amrgenes.fna \
    --plus
    """
}

process amrfinderplus_concat {
    publishDir "${params.outdir}/amrfinderplus", mode: "copy"
    
    label "amrfinderplus"
    
    tag {"Running"}
    
    conda 'bioconda::ncbi-amrfinderplus=3.10.1'
    
    cpus 1

    input:
    path(tsv)
    output:
    path("amrfinderplust_summary.tsv")

    script:
    """
           
    csvtk -t concat ${tsv} > amrfinderplust_summary.tsv
    """
}

process mob_recon {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: { filename -> if (filename.endsWith(".tsv")){ 
                "report/${filename}"  } 
            else if (file(filename).isDirectory())
                {"${filename}"} 
            else { 
                "${filename}"
            }
    }
  
//   container 'kbessonov/mob_suite:3.0.3'

  tag {sample_id}

  conda 'bioconda::mob_suite=3.0.3'
  
  cpus 8
  
  memory '32.GB' 

  input:
  
  tuple val(sample_id), path(contigs)
  
  output:
  path("${sample_id}"), emit: outdir
  tuple val(sample_id), path("${sample_id}/*.fasta"), emit: contigs
  path("*.tsv"), optional: true
  
  script:
  mob_db = params.mob_db ? "-d ${params.mob_db}" : ""
  ref_db = params.ref_db ? "-g ${params.ref_db}" : ""
  """
    hostname > hostname
    export BLASTDB_LMDB_MAP_SIZE=1000000
    mob_recon \
    --num_threads ${task.cpus} \
    --run_typer \
    -s ${sample_id} \
    ${mob_db} \
    ${ref_db} \
    --infile ${contigs} \
    --outdir ${sample_id}
    cp ${sample_id}/contig_report.txt ${sample_id}_contig_report.tsv 2>/dev/null || :
    cp ${sample_id}/mobtyper_results.txt ${sample_id}_mobtyper_results.tsv 2>/dev/null || :
    mv ${sample_id}/chromosome.fasta ${sample_id}/${sample_id}_chromosome.fasta 2>/dev/null || :
  """
}


process CAT {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}

    conda 'bioconda::cat=5.0.0'

    cpus 48

    memory '128.GB'

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}_CAT.*")
    script:
    """
    if [[ \$(file -L ${contigs} | grep -q 'gzip compressed data' && echo yes || echo no) == "yes" ]]
    then
        zcat ${contigs} > ${sample_id}
        
    else
        ln -s ${contigs} ${sample_id}
    fi
    
    CAT contigs -c ${sample_id} -d $params.CAT_db -t $params.CAT_db -n $task.cpus --out_prefix ${sample_id}_CAT

    CAT add_names \
    -i ${sample_id}_CAT.contig2classification.txt \
    -o ${sample_id}_CAT.contig2classification.official_names.txt \
    -t $params.CAT_db \
    --only_official
    """
}


process ABRICATE_RUN {
    label 'abricate'

    publishDir "${params.outdir}/abricate/${db}", mode: 'copy'

    tag {sample_id}
    
    conda 'bioconda::abricate=0.9.7'

    cpus 8

    input:
    tuple val(sample_id), path(contigs), val(db)

    output:
    tuple val(db), path("${sample_id}_${db}.tab")
    
    script:
    """
    abricate --db $db \
    --threads $task.cpus \
    $contigs > ${sample_id}_${db}.tab
    """
}

process ABRICATE_SUMMARY {
    publishDir "${params.outdir}/abricate", mode: 'move'

    label 'abricate'
    
    conda 'bioconda::abricate=0.9.7'

    tag {db}
    
    input:
    tuple val(db), path(tab)

    output:
    path("${db}.tsv")
    
    script:
    """
    abricate --summary $tab > ${db}.tsv
    """
}