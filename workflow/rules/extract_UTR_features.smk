"""
Extract UTR features from different annotations
"""

interval_out_dir = output_root / 'annotations'


rule extract_Gencode_UTR:
    # Extract annotated 3'UTR regions from GENCODE gene annotation
    input: rules.download_Gencode.output  # config["databases"]["Gencode"]["file"]
    output:
        utr=interval_out_dir / 'Gencode/3UTR.bed',
    shell:
        """
        zcat {input} | grep three_prime_UTR | sortBed |\
         mergeBed -c 3,6,7 -o distinct,distinct,distinct > {output.utr}
        # TODO: determine 3UTR start depending on strand
        """


rule extract_PolyA_DB:
    # Extract PAS positions, surrounding regions and hexamers
    input:
        db=rules.download_PolyA_DB.output,# config["databases"]["PolyA_DB"]["file"],
        fasta=ancient(
            expand(rules.genomepy.output[0],assembly=config['assembly_ucsc'])
        )
    output:
        intervals=interval_out_dir / 'PolyA_DB/annotations.tsv',
        stats=interval_out_dir / 'PolyA_DB/stats.tsv'
    params:
        annotations=config['PolyA_DB']['annotation_columns']
    script: '../scripts/extract_polyadb.py'


rule merge_UTR_intervals:
    input:
        utrs=rules.extract_Gencode_UTR.output.utr,  # TODO: use PolyADB UTR intervals
        hexamers=rules.extract_PolyA_DB.output.PAS_hexamers,
        pas=rules.extract_PolyA_DB.output.PAS
    output:
        intervals=interval_out_dir / 'merged_UTR_intervals.tsv'
    params:
        chr_style_gnomAD='' if config['genome_assembly'] == 'GRCh37' else 'chr',
        annotations=INTERVAL_ANNOTATIONS
    script: '../scripts/merge_utr_intervals.py'
