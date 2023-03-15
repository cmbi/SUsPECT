#!/usr/bin/env nextflow

/* 
 * Predict protein function using PolyPhen
 */

process pph2 {
  /*
  Run PolyPhen-2 on a protein sequence with a substitions file

  Returns
  -------
  Returns 2 files:
      1) Output '*.txt'
      2) Errors '*.err'
  */

  tag "${peptide_id}"
  label 'pph2'
  containerOptions "--bind ${params.polyphen2_data}:/opt/pph2/data"
  errorStrategy 'ignore'

  input:
    tuple val(peptide_id), path(fasta), path(subs)

  output:
    tuple val(peptide_id), path("pph2.txt"), path(subs), optional: true

  """
  awk 'BEGIN { OFS="\t" }; {print \$5, \$6, \$7, \$8}' ${subs} > var.subs

  mkdir -p tmp/lock
  run_pph.pl -A -d tmp -s ${fasta} var.subs > pph2.txt

  # Remove output if only contains header
  if [ "\$( wc -l <pph2.txt )" -eq 1 ]; then rm pph2.txt; fi
  """
}

process weka {
  /*
  Run Weka

  Returns
  -------
  Returns 2 files:
      1) Output '*.txt'
      2) Error '*.err'
  */
  tag "${peptide_id} ${model}"
  label 'pph2'
  errorStrategy 'ignore'

  input:
    val model
    tuple val(peptide_id), path(in), path(subs)

  output:
    path 'weka.out'

  shell:
  '''
  run_weka.pl -l /opt/pph2/models/!{model} !{in} 1> weka.txt 2> weka.err

  # clean results and prepend variant coordinates and transcript to each line
  echo "#Chrom,Pos,Ref,Alt,Transcript,Model,Prediction,Score" |\
    sed 's/,/\t/g' > weka.out

  grep -v "^#" weka.txt |\
    paste <(cut -f 1-5 !{subs}) - |\
    awk -v model="!{model}" -v OFS=',' -F'\t' '{gsub(/[ ]{2,}/, "", $0); print $1,$2,$3,$4,$5,model,$17,$21}' |\
    sed 's/,/\t/g' >> weka.out
  '''
}
