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
  label 'vep'

  input:
    file weka_out
    file sift_out
    path vep_config
    path dir_plugins

  output:
    path 'vep_config_plugin.ini'

  """
  new_config="vep_config_plugin.ini"
  cp ${vep_config} \${new_config}
  echo "dir_plugins \$(realpath ${dir_plugins})" >> \${new_config}

  # Prepare PolyPhen-2 results
  weka=\$(realpath ${weka_out})
  sort -k1 -nk2 \${weka} | uniq | bgzip > \${weka}.gz
  tabix \${weka}.gz -b 2 -e 2
  echo "plugin TranscriptAnnotator,file=\${weka}.gz,prefix=PolyPhen2_" >> \${new_config}

  # Prepare SIFT results
  sift=\$(realpath ${sift_out})
  sort -k1 -nk2 \${sift} | uniq | bgzip > \${sift}.gz
  tabix \${sift}.gz -b 2 -e 2
  echo "plugin TranscriptAnnotator,file=\${sift}.gz,prefix=SIFT_" >> \${new_config}
  """
}

process create_exclusion_variants {
  label 'vep'
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
  label 'bedtools'
  storeDir "${params.outdir}/${params.name}/filtering/"

  input:
    file vcf
    file exclusion_vcf

  output:
    path 'benign.vcf'

  """
  bedtools intersect -v -a $vcf -b $exclusion_vcf -wa -header > benign.vcf
  """
}

process filter_common_variants {
  label 'vep'
  storeDir "${params.outdir}/${params.name}/filtering/"

  input:
    file vcf

  output:
    path 'vep_filtered.vcf', emit: vep_filtered_vcf
    path 'benign.vcf.gz', emit: originally_benign_af_vcf
    path 'benign.vcf.gz.tbi', emit: originally_benign_af_vcf_index

  """
  filter_vep -i $vcf --filter "MAX_AF < 0.01 or not MAX_AF" | filter_vep -o vep_filtered.vcf --only_matched --filter "Feature matches _ORF_" --filter "Consequence is missense_variant"
  bgzip $vcf
  tabix -p vcf ${vcf}.gz
  """
}

// process filter_common_variants {
//   label 'vep'

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
