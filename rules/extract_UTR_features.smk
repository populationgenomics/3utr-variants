"""
Extract UTR features from different annotations
"""

interval_out_dir = output_root / 'intervals'

rule extract_Gencode_UTR:
    input: rules.download_Gencode.output  # config["databases"]["Gencode"]["file"]
    output:
        utr=interval_out_dir/'Gencode/3UTR.bed',
        pas=interval_out_dir/'Gencode/PAS.bed'
    shell:
        """
        zcat {input} | grep three_prime_UTR | gff2bed > {output.utr}
        # TODO: determine 3UTR start depending on strand
        """


rule extract_PolyA_DB:
    input:
        db=rules.download_PolyA_DB.output,  # config["databases"]["PolyA_DB"]["file"],
        fasta=ancient(
            expand(rules.genomepy.output[0], assembly=config['assembly_ucsc'])
        )
    output:
        PAS=interval_out_dir/'PolyA_DB/PAS-{filter}.bed',
        PAS_context_40nt=interval_out_dir/'PolyA_DB/context_40nt-{filter}.bed',
        PAS_context_100nt=interval_out_dir/'PolyA_DB/context_100nt-{filter}.bed',
        PAS_hexamers=interval_out_dir/'PolyA_DB/hexamers-{filter}.bed',
        stats=interval_out_dir/'PolyA_DB/stats-{filter}.txt'
    script: '../scripts/extract_polyadb.py'
