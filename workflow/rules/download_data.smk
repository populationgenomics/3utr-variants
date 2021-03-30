"""
Download genome sequence reference, gene and 3'UTR annotation files
"""


rule genomepy:
    # Download reference sequence using genomepy
    output:
        protected(
            multiext(str(output_root / "reference/{assembly}/{assembly}"),
                ".fa",".fa.fai",".fa.sizes",".annotation.gtf.gz",".annotation.bed.gz")
        )
    log: str(output_root / "logs/genomepy_{assembly}.log")
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
    # Request fasta file for genome assembly specified in config
    input: expand(rules.genomepy.output,assembly=config['assembly_ucsc'])


rule download_PolyA_DB:
    # Download and extract PolyA_DB v3 annotation
    output: output_root / 'annotations/PolyA_DB/human.PAS.txt'
    params:
        url=config['PolyA_DB']['url']
    shell:
        """
        tmpdir="$(dirname "$(tempfile)")"
        # TODO: download from GCP bucket
        wget --no-check-certificate -nc -P $tmpdir {params.url}
        unzip "$tmpdir"/human_pas.zip -d $(dirname {output})
        """


rule download_Gencode:
    # Download GENCODE gene annotation
    output: output_root / 'annotations/Gencode/gencode.v36lift37.annotation.gff3.gz'
    params:
        url=config['Gencode']['url']
    shell:
        """
        wget -nc -P $(dirname {output}) {params.url}
        """


rule download_PolyASite2:
    output: output_root / 'annotations/PolyASite2/atlas.clusters.2.0.GRCh38.96.bed'
    params:
        url = config['PolyASite2']['url']
    shell:
        """
        wget -nc -P $(dirname {output}) {params.url}
        gzip -d {output}.gz
        """
