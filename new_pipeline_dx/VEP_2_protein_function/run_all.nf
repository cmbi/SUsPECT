/*
 * From VCF to protein prediction
 */

nextflow.enable.dsl=2

params.help = false
params.outdir = "outdir"

// print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP and SIFT/PolyPhen2 based on a VCF file'
  log.info '----------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow run run_all.nf --fasta $fasta --gtf $gtf --vcf $vcf --chros 21,22 -profile lsf -resume'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf VCF      VCF containing missense variants for protein function prediction'
  exit 1
}

include { run_vep } from '../VEP/workflows/run_vep.nf'
include { getTranslation } from '../protein_function/nf_modules/run_agat.nf'
include { pph2 } from '../protein_function/nf_modules/run_polyphen2.nf'
include { create_subs } from './nf_modules/create_subs.nf'
include { linearise_fasta; get_fasta } from './nf_modules/fasta.nf'

workflow {
  // Get translated FASTA
  getTranslation( file(params.gtf), file(params.fasta) )
  linearise_fasta( getTranslation.out.collect() )

  // Get substitutions
  run_vep()  
  create_subs( run_vep.out.vcfFile )
  subs = create_subs.out.flatten()
  translated = get_fasta ( linearise_fasta.out, subs )

  // For each transcript, predict protein function
  pph2(translated.fasta, subs)
}
