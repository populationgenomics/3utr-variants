"""
Extract UTR features from different annotations
"""

interval_out_dir = output_root / 'annotations'


def chr_style_gnomAD(wildcards):
    return '' if config['genome_assembly'] == 'GRCh37' else 'chr'


rule extract_Gencode_UTR:
    # Extract annotated 3'UTR regions from GENCODE gene annotation
    input: rules.download_Gencode.output
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
        db=rules.download_PolyA_DB.output,
        fasta=ancient(
            expand(rules.genomepy.output[0],assembly=config['assembly_ucsc'])
        )
    output:
        intervals=interval_out_dir / 'PolyA_DB/annotations.tsv',
        stats=interval_out_dir / 'PolyA_DB/stats.tsv'
    params:
        annotations=config['PolyA_DB']['annotation_columns']
    script: '../scripts/extract_polyadb.py'


rule extract_PolyASite2:
    input:
        db=rules.download_PolyASite2.output,
        genes=expand(rules.genomepy.output[4],assembly='hg38')
    output:
        intervals=interval_out_dir / 'PolyASite2/annotation.tsv',
        stats=interval_out_dir / 'PolyASite2/stats.tsv'
    params:
        annotations=config['PolyASite2']['annotation_columns']
    script: '../scripts/extract_polyasite2.py'


rule merge_UTR_intervals:
    # Merge all interval annotations with 3'UTR regions
    input:
        utrs=rules.extract_Gencode_UTR.output.utr,
        PolyA_DB=rules.extract_PolyA_DB.output.intervals,
        PolyASite2=rules.extract_PolyASite2.output.intervals,
        chainfile=rules.download_chainfile.output
    output:
        intervals=interval_out_dir / 'merged_UTR_intervals.tsv',
        bed=interval_out_dir / 'merged_UTR_intervals.bed'
    params:
        chr_style_gnomAD=chr_style_gnomAD,
        annotations=ALL_ANNOTATIONS
    script: '../scripts/merge_utr_intervals.py'
