process filter_first_vep {
    publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
    cpus 1
    label 'vcfparser'

    input:
      path gff3_file
      path vep_output
    output:
      path "filtered_file.vcf"
      path "Variants_more_severe_in_new_annotation.csv"
      

    script:
        """
        sed -i '1d' $gff3_file
        python $workflow.projectDir/bin/filter_vep.py $gff3_file $vep_file filtered_file.vcf
        """

}
