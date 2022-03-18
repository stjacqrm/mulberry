#!/usr/bin/env nextflow

//Description: Mulberry is a pipeline designed to analyse and visualize AMR gene frequency data
//Author: Rachael St. Jacques
//email: rachael.stjacques@dgs.virginia.gov


//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
    .set { raw_reads }

Channel
    .fromPath("$baseDir/report/report_template.Rmd")
    .set { report }

Channel
    .fromPath("$baseDir/bash/create_mulberry_report.sh")
    .set { bash }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(reads) into read_files_fastqc, read_files_trimming

  script:
  if(params.name_split_on!=""){
    name = name.split(params.name_split_on)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    """
    """
      }
  }

//Step1a: FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/logs/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from read_files_fastqc

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  script:
  """
  fastqc -q  ${reads}
  """
  }

//Step1b: Trim with Trimmomatic
process trim {
  tag "$name"
  if(params.savetrimmedreads){
    publishDir "${params.outdir}/trimmed", mode: 'copy'
  }
  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_trimmed{_1,_2}.fastq.gz") into trimmed_reads
  file("${name}.trim.stats.txt") into trimmomatic_stats

  script:
  """
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads ${task.cpus} ${reads} -baseout ${name}.fastq.gz SLIDINGWINDOW:${params.windowsize}:${params.qualitytrimscore} MINLEN:${params.minlength} 2> ${name}.trim.stats.txt
  mv ${name}*1P.fastq.gz ${name}_trimmed_1.fastq.gz
  mv ${name}*2P.fastq.gz ${name}_trimmed_2.fastq.gz
  """
}

//Step2: Remove PhiX contamination
process cleanreads {
  tag "$name"
  publishDir "${params.outdir}/logs/cleanedreads/stats", mode: 'copy',pattern:"*.stats.txt"
  publishDir "${params.outdir}/logs/cleanedreads/reads", mode: 'copy',pattern:"*.fastq.gz"

  input:
  set val(name), file(reads) from trimmed_reads

  output:
  tuple name, file("${name}{_1,_2}.clean.fastq.gz") into cleaned_reads, reads_files_minimap
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats

  script:
  """
  repair.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh -Xmx"${task.memory.toGiga()}g" in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_1.clean.fastq.gz out2=${name}_2.clean.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  """
}


//Step3: Assemble trimmed reads with Shovill
process shovill {
  errorStrategy 'ignore'
  tag "$name"
  publishDir "${params.outdir}/results/assemblies", mode: 'copy'

  input:
  set val(name), file(reads) from cleaned_reads

  output:
  tuple name, file("${name}.fasta") into assembled_genomes_quality, assembled_genomes,assembled_genomes_ar,assembled_genomes_gamma

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory}  --outdir . --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv contigs.fa ${name}.fasta
  """
}

//Step3a: Assembly Quality Report
process quast {
  errorStrategy 'ignore'
  publishDir "${params.outdir}/logs/quast",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_report

  script:
  """
  quast.py ${assembly} -o .
  mv report.txt ${name}.quast.tsv
  """
}


//Step4a: Find AR genes with amrfinder+
process amrfinder {
  tag "$name"
  publishDir "${params.outdir}/results/amrfinder",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_ar

  output:
  file("${name}.tsv") into ar_predictions

  script:
  """
  amrfinder -n ${assembly} -o ${name}.tsv
  """
}

//temp Step4b: Find AR genes with amrfinder+
process concat_amr_files {
  //tag "$name"
  publishDir "${params.outdir}/results",mode:'copy'

  input:
  file(ar) from ar_predictions.collect()

  output:
  file("ar_prediction.tsv") into ar_output

  script:
  """
  awk '(NR == 1) || (FNR > 1)' *.tsv > ar_prediction.tsv
  """
}
//Step4b: Find AR genes with GAMMA
/*process gamma {
  tag "$name"
  publishDir "${params.outdir}/results/gamma",mode:'copy'

  input:
  set val(name), file(assembly) from assembled_genomes_gamma

  output:
  file("${name}.tsv") into gamma_predictions

  script:
  """
  for i in *; do GAMMA.py "\$i" ~/Applications/GAMMA/Five_Vir_markers_all.fasta "\$i"_gamma; mkdir gamma; mv "\$i"_gamma.gamma "\$i"_gamma.tsv; mv *tsv gamma; mkdir psl; mv *.psl psl; done

  """
}*/



//Step 6: generate the pdf report
process render{
  publishDir "${params.outdir}/", mode: 'copy', pattern:'Mulberry-report.pdf'
  publishDir "${params.outdir}/", mode: 'copy' , pattern:"*csv"
  beforeScript 'ulimit -s unlimited'

  input:
  file(ar) from ar_output
  file(rmd) from report
  file(bash) from bash

  output:
  file('*.csv')
  file "Mulberry-report.pdf"

  shell:
  """
  cp ${rmd} ./report_template.Rnw
  cp ${bash} ./create_report.sh
  chmod +x create_report.sh
  bash create_report.sh -p "Mulberry-report" -t "${params.title}" -T report_template.Rnw -o . -a ${ar}
  """

}
