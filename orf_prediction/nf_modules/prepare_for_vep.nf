process gtf_for_vep{
    storeDir "${params.outdir}/${params.name}/final_gtf/"
    cpus 1

    label 'vep'

    input:
        path complete_gtf
    output:
        path "novel_complete.gtf.gz"
        path "novel_complete.gtf.gz.tbi"

    script:
        """
        sed 's/; \$/;/g' $complete_gtf | grep -v "#" | sort -k1,1 -k4,4n -k5,5n -t\$'\t' | bgzip -c > novel_complete.gtf.gz
        tabix -p gff novel_complete.gtf.gz
        """
}
