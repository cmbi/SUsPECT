process gtf_for_vep{
    publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
    cpus 1
    container 'quay.io/biocontainers/ensembl-vep:106.1--pl5321h4a94de4_0'

    input:
        path complete_gff
    output:
        path "novel_complete.gff.gz"
        path "novel_complete.gff.gz.tbi"

    script:
        """
        grep -v "#" $complete_gff | sort -k1,1 -k4,4n -k5,5n -t\$'\t' | bgzip -c > novel_complete.gff.gz
        tabix -p gff novel_complete.gff.gz
        """
}