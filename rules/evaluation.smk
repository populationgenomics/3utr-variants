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


def maps_files(wildcards):
    subsets = variant_subsets.keys()
    if wildcards.run_location == 'local':
        return expand(rules.MAPS_local.output.maps, variant_subset=subsets)
    return GS.remote(expand(rules.MAPS_GCP.output.maps.__str__(), variant_subset=subsets))


rule gather_MAPS:
    input: maps_files
    output: output_root / 'MAPS/all_{run_location}.tsv'
    run:
        import pandas as pd

        df_list = []
        for variant_subset, file in zip(variant_subsets.keys(), input):
            df = pd.read_table(file.__str__(), sep='\t')
            df['variant_subset'] = variant_subset
            df_list.append(df)
        df_all = pd.concat(df_list)
        df_all.to_csv(output[0], index=False, sep='\t')

rule plots:
    input: rules.gather_MAPS.output
    output:
          maps=output_root / 'plots/MAPS_{run_location}.png'
    script: '../scripts/plots.py'
