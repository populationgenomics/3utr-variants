"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()


variant_subsets = {
    'PolyADB_40nt-conserved': expand(
        rules.extract_PolyA_DB.output.PAS_context_40nt,
        filter='conserved'
    ),
    'PolyADB_100nt-conserved': expand(
        rules.extract_PolyA_DB.output.PAS_context_100nt,
        filter='conserved'
    ),
    #'PolyADB_hexamers-conserved': expand(
    #    rules.extract_PolyA_DB.output.PAS_hexamers,
    #    filter='conserved'
    #),
    'PolyADB_hexamers-all': expand(
        rules.extract_PolyA_DB.output.PAS_hexamers,
        filter='all'
    ),
    'Gencode_UTR': rules.extract_Gencode_UTR.output.utr
}

rule convert_bedfile:
    input: lambda wildcards: variant_subsets[wildcards.variant_subset]
    output: output_root/'intervals/{variant_subset}_hail.txt'
    params:
        chr_style_hail = '' if config['genome_assembly'] == 'GRCh37' else 'chr'
    script: '../scripts/convert_bed.py'

rule MAPS_GCP:
    input: rules.convert_bedfile.output
    output:
        maps=GS.remote(
            f'{config["bucket"]}/MAPS_{{variant_subset}}.tsv',
            keep_local=True
        ),
    log: GS.remote(f'{config["bucket"]}/MAPS_{{variant_subset}}.log')
    params:
        gnomad_prepare='scripts/prepare_gnomad.py',
        maps_score='scripts/maps_score.py'
    shell:
        """
        hailctl dataproc submit {config[cluster]} \
            --files {input} \
            --pyfiles {params.gnomad_prepare},{params.maps_score} \
            scripts/maps_gcp.py  \
                -o gs://{output.maps} --log gs://{log} \
                --intervals $(basename {input}) \
                --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
                --context_ht {config[gnomAD][context_ht]} \
                --mutation_ht {config[gnomAD][mutation_rate_ht]} \
                --genome_assembly {config[genome_assembly]}
        """

rule prepare_gnomAD:
    threads: 10
    output:
        gnomAD_ht=directory(output_root/'gnomAD.ht')
    log: hail=str(output_root/'logs/prepare_gnomAD_hail.log')
    script: '../scripts/prepare_gnomad.py'

rule MAPS_local:
    input:
         intervals=rules.convert_bedfile.output,
         gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
        maps=output_root/'MAPS/{variant_subset}_local.tsv',
    log: hail=str(output_root/'logs/MAPS_local_{variant_subset}_hail.log')
    threads: 10
    script: '../scripts/maps_score.py'

rule plots:
    input:
        expand(
            rules.MAPS_local.output if config['local'] else rules.MAPS_GCP.output,
            variant_subset=variant_subsets.keys()
        )
    output:
        maps=output_root/'plots/MAPS.png'
    script: '../scripts/plots.py'
