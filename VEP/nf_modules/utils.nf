process append_fasta_gtf_to_config {
  input:
    path vep_config
    path fasta
    path gtf

  output:
    path '*.ini'

  """
  new_config="vep_config_complete.ini"
  cp ${vep_config} \${new_config}

  echo "max_af 1" >> \${new_config}
  echo "polyphen p" >> \${new_config}
  echo "gtf \$(realpath ${gtf})" >> \${new_config}
  echo "fasta \$(realpath ${fasta})" >> \${new_config}
  """
}

process prepare_vep_transcript_annotation {
  container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

  input:
    file weka_out
    path vep_config
    path dir_plugins

  output:
    //path '*.ini'
    path '*_plugin.ini'

  """
  new_config="vep_config_plugin.ini"
  cp ${vep_config} \${new_config}

  weka=\$(realpath ${weka_out})
  sort -V -k1,1 -k2,2 \${weka} | bgzip > \${weka}.gz
  tabix \${weka}.gz -b 2 -e 2

  echo "plugin TranscriptAnnotator,\${weka}.gz" >> \${new_config}
  echo "dir_plugins \$(realpath ${dir_plugins})" >> \${new_config}
  """
}

process create_exclusion_variants {
  container 'ensemblorg/ensembl-vep:latest'
  storeDir "${params.outdir}/${params.name}/filtering/"

  input:
    file vcf

  output:
    path 'exclude.vcf'

  """
  zgrep "^#" $vcf > vep_orf.vcf 
  zgrep -v "^#" $vcf | grep "_ORF_" >> vep_orf.vcf
  filter_vep -i vep_orf.vcf --only_matched --filter "Feature not matches _ORF_" | filter_vep -o exclude.vcf --filter "PolyPhen matches damaging or IMPACT is HIGH"
  """
}

process exclude_pathogenic {
  container 'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0'
  storeDir "${params.outdir}/${params.name}/filtering/"

  input:
    file vcf
    file exclusion_vcf

  output:
    path 'benign.vcf'

  """
  bedtools intersect -v -a $vcf -b $exclusion_vcf -wa > benign.vcf
  """
}

process filter_common_variants {
  container 'ensemblorg/ensembl-vep:latest'
  storeDir "${params.outdir}/${params.name}/filtering/"

  input:
    file vcf
    file old_vcf

  output:
    path 'vep_filtered.vcf', emit: vep_filtered_vcf
    path 'originally_benign_af.vcf.gz', emit: originally_benign_af_vcf
    path 'originally_benign_af.vcf.gz.tbi', emit: originally_benign_af_vcf_index

  """
  zgrep "^#" $old_vcf > originally_benign.vcf
  cat $vcf >> originally_benign.vcf
  filter_vep -i originally_benign.vcf -o originally_benign_af.vcf --filter "MAX_AF < 0.01 or not MAX_AF" 
  filter_vep -i originally_benign_af.vcf -o vep_filtered.vcf --only_matched --filter "Feature matches _ORF_" --filter "Consequence is missense_variant"
  bgzip originally_benign_af.vcf
  tabix -p vcf originally_benign_af.vcf.gz
  """
}

// process filter_common_variants {
//   container 'ensemblorg/ensembl-vep:latest'

//   input:
//     file vcf

//   output:
//     path 'vep_filtered.vcf'

//   """
//   zgrep "^#" $vcf > vep_orf.vcf 
//   zgrep -v "^#" $vcf | grep -e "_ORF_" -e "missense" >> vep_orf.vcf
//   filter_vep -i vep_orf.vcf --only_matched --filter "Feature matches _ORF_" --filter "Consequence is missense_variant" | filter_vep --filter "MAX_AF < 0.01 or not MAX_AF" -o vep_filtered.vcf
//   """
// }