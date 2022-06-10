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

  echo "gtf \$(realpath ${gtf})" >> \${new_config}
  echo "fasta \$(realpath ${fasta})" >> \${new_config}
  echo "max_af" >> \${new_config}
  """
}

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

process filter_common_variants {
  container 'ensemblorg/ensembl-vep:latest'

  input:
    file vcf

  output:
    path 'vep_filtered.vcf'

  """
  grep "^#" > vep_orf.vcf 
  grep -v "^#" | grep "_ORF_" >> vep_orf.vcf
  filter_vep -i vep_orf.vcf -o vep_filtered.vcf --filter "AF < 0.01 or not AF"
  """
}
