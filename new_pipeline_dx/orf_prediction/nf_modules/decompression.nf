/*--------------------------------------------------
Decompress Logit Model
---------------------------------------------------*/
process gunzip_logit_model {
   tag "decompress logit model"
   cpus 1

   input:
   path logit_model

   output:
   path "*.RData"

   script:
   """
   gunzip -f ${logit_model}
   """
}

/*--------------------------------------------------
Decompress VCF file
---------------------------------------------------*/
process gunzip_vcf {
   tag "decompress VCF file"
   cpus 1

   input:
   path vcf

   output:
   path "*.vcf"

   script:
   """
   gunzip -f ${vcf}
   """
}

/*--------------------------------------------------
Decompress genome fasta file
---------------------------------------------------*/
process gunzip_genome_fasta {
   tag "decompress gzipped genome fasta"
   cpus 1

   input:
   path genome_fasta

   output:
   path "*.{fa,fasta}"

   script:
   """
   gunzip -f ${genome_fasta}
   """
}