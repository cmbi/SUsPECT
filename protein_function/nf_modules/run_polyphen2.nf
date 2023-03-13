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

  tag "${subs.chrom}-${subs.pos}-${subs.ref}-${subs.alt}-${subs.feature}"
  label 'pph2'
  containerOptions "--bind ${params.polyphen2_data}:/opt/pph2/data"
  errorStrategy 'ignore'

  input:
    path fasta
    val subs

  output:
    path "*.txt", emit: results

  """
  echo "${subs.feature}\t${subs.protein_pos}\t${subs.aa_ref}\t${subs.aa_alt}" > var.subs
  grep -A1 ${subs.feature} ${fasta} > peptide.fa

  mkdir -p tmp/lock
  out="${subs.chrom}-${subs.pos}-${subs.ref}-${subs.alt}-${subs.feature}.txt"
  run_pph.pl -A -d tmp -s peptide.fa var.subs > \$out

  # Remove output if only contains header
  if [ "\$( wc -l <\$out )" -eq 1 ]; then rm \$out; fi
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
  tag "${in.baseName}"
  label 'pph2'
  errorStrategy 'ignore'

  input:
    val model
    path in

  output:
    path '*.txt', emit: results
    path '*.out', emit: processed
    path '*.err', emit: error

  shell:
  '''
  res=!{in.baseName}_!{model}.txt
  run_weka.pl -l /opt/pph2/models/!{model} !{in} \
              1> $res 2> !{in.baseName}_!{model}.err

  # clean results and append variant coordinates and transcript id
  var=$( echo !{in.baseName}-PolyPhen2 | sed 's/-/,/g' )
  grep -v "^#" ${res} | awk -v var="$var" -v OFS=',' -F'\t' '{$1=$1;print var,$12,$16}' | sed 's/,/\t/g' > !{in.baseName}.out
  '''
}
