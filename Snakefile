configfile: "config.yml"
output_root = config["output_root"]

rule all:
    input: output_root + "/.done"

rule dependency:
    output: output_root + "/dependency.svg"
    shell: "snakemake --dag | dot -Tsvg -Grankdir=TB > {output}"

rule prepare_gnomAD:
    threads: 10
    output:
        checkpoint=temp(directory(output_root + "/gnomAD_tmp.ht")),
        gnomAD_ht=directory(output_root + "/gnomAD.ht")
    script: "scripts/prepare_gnomAD.py"

rule extract_Gencode_UTR:
    input: config["databases"]["Gencode"]["file"]
    output: output_root + "/Gencode/3UTR.bed"
    shell:
        """
        zcat {input} | grep three_prime_UTR | gff2bed > {output}
        """

rule extract_Gencode_PAS:
    input: rules.extract_Gencode_UTR.output
    output: output_root + "/Gencode/PAS.gff3"
    run:
        """
        # TODO: determine 3UTR start depending on strand
        """

rule extract_PolyA_DB:
    input: config["databases"]["PolyA_DB"]["file"]
    output:
        PAS=output_root + "/PolyA_DB/PAS.bed",
        PAS_context_40nt=output_root + "/PolyA_DB/PAS_context_40nt.bed",
        PAS_context_100nt=output_root + "/PolyA_DB/PAS_context_100nt.bed"
    script: 'scripts/extract_PolyA_DB.py'


rule score_MAPS:
    # TODO: generalise by DB type
    input:
         bed=rules.extract_Gencode_UTR.output,
         gnomAD=rules.prepare_gnomAD.output['gnomAD_ht']
    output:
          maps=directory(output_root + "/Gencode/MAPS.ht"),
          done=touch(output_root + "/.done")
    threads: 10
    script: "scripts/maps_score.py"
