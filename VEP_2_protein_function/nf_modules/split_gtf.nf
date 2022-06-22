process split_gtf {
  /*
  Split GTF into multiple files
  */

  tag "${gtf.baseName}"
  container "quay.io/biocontainers/agat:0.9.0--pl5321hdfd78af_0"

  input:
    path gtf

  output:
    path '*.gtf'

  """
  agat_sq_split.pl --input $gtf -o split.gtf
  """
}
