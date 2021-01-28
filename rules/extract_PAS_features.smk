output_root = config["output_root"]

rule extract_Gencode_UTR:
    input: config["databases"]["Gencode"]["file"]
    output: output_root + "/Gencode/3UTR.bed"
    shell:
        """
        zcat {input} | grep three_prime_UTR | gff2bed > {output}
        # TODO: manage chromosome style depending on gene assembly
        """

rule extract_Gencode_PAS:
    input: rules.extract_Gencode_UTR.output
    output: output_root + "/Gencode/PAS.gff3"
    run:
        """
        # TODO: determine 3UTR start depending on strand
        """

rule extract_PolyA_DB:
    input:
        db=config["databases"]["PolyA_DB"]["file"],
        fasta=expand(rules.genomepy.output[0], assembly=config['assembly_ucsc'])
    output:
        PAS=output_root + "/PolyA_DB/PAS.bed",
        PAS_context_40nt=output_root + "/PolyA_DB/PAS_context_40nt.bed",
        PAS_context_100nt=output_root + "/PolyA_DB/PAS_context_100nt.bed",
        PAS_hexamers=output_root + "/PolyA_DB/PAS_hexamers.bed"
    script: '../scripts/extract_PolyA_DB.py'
