output_root = config["output_root"]

rule genomepy:
    output:
        multiext(output_root + "/reference/{assembly}/{assembly}",
                 ".fa", ".fa.fai", ".fa.sizes", ".annotation.gtf.gz", ".annotation.bed.gz")
    log:
        output_root + "/logs/genomepy_{assembly}.log"
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

assembly_ucsc = 'hg19' if config['genome_assembly'] == 'GChr37' else 'GChr38'

rule get_fasta:
    input: expand(rules.genomepy.output, assembly=assembly_ucsc)
