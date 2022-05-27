process gtf_for_vep{
    publishDir "${params.outdir}/${params.name}/final_gtf/", mode: 'copy'
    cpus 1
    conda 'bioconda::pbgzip bioconda::tabix'

    input:
        path complete_gtf
    output:
        path "${params.name}_complete.gtf.gz"

    script:
        """
        grep -v "#" $complete_gtf | sort -k1,1 -k4,4n -k5,5n -t\$'\t' | bgzip -c > ${complete_gtf}.gz
        tabix -p gff ${complete_gtf}.gz
        """
}