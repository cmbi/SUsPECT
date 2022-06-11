#!/usr/bin/env nextflow

/* 
 * Predict protein function using PolyPhen
 */

params.polyphen2_data = "/hps/nobackup/flicek/ensembl/variation/nuno/sift-polyphen2-nextflow-4667/input/polyphen2"

// PolyPhen-2 data directories to bind
PPH="/opt/pph2" // PolyPhen-2 directory in Singularity container
dssp="${params.polyphen2_data}/dssp:$PPH/dssp"
wwpdb="${params.polyphen2_data}/wwpdb:$PPH/wwpdb"
precomputed="${params.polyphen2_data}/precomputed:$PPH/precomputed"
nrdb="${params.polyphen2_data}/nrdb:$PPH/nrdb"
pdb2fasta="${params.polyphen2_data}/pdb2fasta:$PPH/pdb2fasta"
ucsc="${params.polyphen2_data}/ucsc:$PPH/ucsc"
uniprot="${params.polyphen2_data}/uniprot:$PPH/uniprot"

process pph2 {
  /*
  Run PolyPhen-2 on a protein sequence with a substitions file

  Returns
  -------
  Returns 2 files:
      1) Output '*.txt'
      2) Errors '*.err'
  */

  tag "${subs.baseName}"
  container "nunoagostinho/polyphen-2:2.2.3"
  containerOptions "--bind $dssp,$wwpdb,$precomputed,$nrdb,$pdb2fasta,$ucsc,$uniprot"
  memory '4 GB'
  errorStrategy 'ignore'

  input:
    path protein
    path subs

  output:
    path '*.txt', emit: results
    path '*.err', emit: error

  publishDir "${params.outdir}/polyphen-2"

  """
  mkdir -p tmp/lock
  output=${subs.baseName}.txt
  run_pph.pl -A -d tmp -s $protein $subs \
             1> \$output 2> ${subs.baseName}.err
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
  container "nunoagostinho/polyphen-2:2.2.3"
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
