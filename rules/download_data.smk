rule genomepy:
    output:
        protected(
            multiext(str(output_root/"reference/{assembly}/{assembly}"),
                 ".fa", ".fa.fai", ".fa.sizes", ".annotation.gtf.gz", ".annotation.bed.gz")
        )
    log: str(output_root/"logs/genomepy_{assembly}.log")
    params:
        provider="UCSC"  # optional, defaults to ucsc. Choose from ucsc, ensembl, and ncbi
    cache: True  # mark as eligible for between workflow caching
    shell:
        """
        genome_dir=$(dirname "{output[0]}")
        genome_dir=$(dirname "$genome_dir")
        echo downloading genome to: $genome_dir
        genomepy install {wildcards.assembly} \
            -p {params.provider} --annotation \
            -g $genome_dir >> {log} 2>&1
        """

rule get_fasta:
    input: expand(rules.genomepy.output, assembly=config['assembly_ucsc'])

# TODO: download PolyA_DB3
