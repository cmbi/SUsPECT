/*
 * From VCF to protein prediction
 */

nextflow.enable.dsl=2

params.help = false
params.outdir = "outdir"
params.gff = null
params.vep_config = "../VEP/nf_config/vep.ini"

// print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to run VEP and SIFT/PolyPhen2 based on a VCF file'
  log.info '----------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow run run_all.nf --fasta $fasta --gff $gff --vcf $vcf --chros 21,22 -profile lsf -resume'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf VCF      VCF containing missense variants for protein function prediction'
  exit 1
}

include { run_vep; run_vep as run_vep_plugin } from '../VEP/workflows/run_vep.nf'
include { getTranslation } from '../protein_function/nf_modules/run_agat.nf'
include { pph2; weka } from '../protein_function/nf_modules/run_polyphen2.nf'
include { create_subs } from './nf_modules/create_subs.nf'
include { linearise_fasta; get_fasta } from './nf_modules/fasta.nf'

process prepare_vep_transcript_annotation {
  container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

  input:
    file weka_out
    path vep_config
    path dir_plugins

  output:
    path '*.ini'

  """
  new_config="vep_config_plugin.ini"
  cp ${vep_config} \${new_config}

  weka=\$(realpath ${weka_out})
  sort \${weka} | bgzip > \${weka}.gz
  tabix \${weka}.gz -b 2 -e 2

  echo "plugin TranscriptAnnotator,\${weka}.gz" >> \${new_config}
  echo "dir_plugins \$(realpath ${dir_plugins})" >> \${new_config}
  """
}

workflow {
  // Get translated FASTA
  getTranslation( file(params.gff), file(params.fasta) )
  linearise_fasta( getTranslation.out.collect() )

  // Get substitutions
  vep_config = file( params.vep_config )
  run_vep( vep_config )
  create_subs( run_vep.out.vcfFile )
  subs = create_subs.out.flatten()
  translated = get_fasta ( linearise_fasta.out, subs )

  // For each transcript, predict protein function
  pph2(translated.fasta, subs)
  weka("HumDiv.UniRef100.NBd.f11.model", pph2.out.results)
  res = weka.out.processed.collectFile(name: "weka_results.out")

  // Incorporate PolyPhen-2 scores into VEP
  prepare_vep_transcript_annotation( res, vep_config, file("VEP_plugins") )
  run_vep_plugin( prepare_vep_transcript_annotation.out )
}
