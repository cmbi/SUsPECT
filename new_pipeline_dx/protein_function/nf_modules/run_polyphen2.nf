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
    path '*.txt'
    path '*.err'

  publishDir "${params.outdir}/polyphen-2"

  """
  mkdir -p tmp/lock
  run_pph.pl -A -d tmp -s $protein $subs \
             1> ${subs.baseName}.txt 2> ${subs.baseName}.err
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

  container "nunoagostinho/polyphen-2:2.2.3"
  errorStrategy 'ignore'

  input:
    path model
    path in

  output:
    path '${model}.*'

  """
  run_weka.pl -l $model $in 1> ${model}.txt 2> ${model}.err
  """
}
