#!/usr/bin/env nextflow

/* 
 * Script to merge chromosome-wise VCF files into single VCF file
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
mergedVCF = "merged-file"
params.outdir = ""
params.cpus = 1

process mergeVCF {
  /*
  Function to merge chromosome-wise VCF files into single VCF file

  Returns
  -------
  Returns 2 files:
      1) A VCF format file 
      2) A tabix index for that VCF
  */

  publishDir "${params.outdir}/${params.name}", 
    enabled: "${params.outdir}/${params.name}" as Boolean,
    mode:'copy'
    
  cpus params.cpus
  label 'bcftools'
   
  input:
  path(vcfFiles)
  path(indexFiles)

  output:
  path("${ mergedVCF }.vcf.gz"), emit: vcfFile                                     
  path("${ mergedVCF }.vcf.gz.tbi"), emit: indexFile

  script: 
  """
  bcftools concat ${ vcfFiles } -Oz -o temp-${ mergedVCF}.vcf.gz
  bcftools sort -Oz temp-${ mergedVCF}.vcf.gz -o ${ mergedVCF}.vcf.gz 
  bcftools  index -t ${ mergedVCF}.vcf.gz
  """
}
