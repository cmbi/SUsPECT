process split_gtf {
  /*
  Split GTF into multiple files
  */

  tag "${gtf.baseName}"
  label 'agat'

  input:
    path gtf

  output:
    path '*.gtf'

  """
  agat_sq_split.pl --input $gtf -o split.gtf
  """
}
