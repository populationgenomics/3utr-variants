"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""

rule convert_bedfile:
    input: lambda wildcards: variant_subsets[wildcards.variant_subset]
    output: output_root/'intervals/{variant_subset}_hail.txt'
    params:
        chr_style_hail = '' if config['genome_assembly'] == 'GRCh37' else 'chr'
    script: '../scripts/convert_bed.py'

rule MAPS_GCP:
    input: rules.convert_bedfile.output
    output:
        maps=directory(output_root/'MAPS/{variant_subset}.ht')
    log:
        local=str(output_root/'logs/MAPS_GCP_{variant_subset}_hail.log')
    params:
        gnomad_prepare='scripts/prepare_gnomad.py',
        maps_score='scripts/maps_score.py'
    shell:
        """
        OUTPUT_GCP="maps_{wildcards.variant_subset}.ht"
        LOG_GCP="gs://{config[bucket]}/MAPS_GCP_{wildcards.variant_subset}.log"
        # TODO: manage GCP files with snakemake.remote.GS remote provider
        hailctl dataproc submit {config[cluster]} \
            --files {input} \
            --pyfiles {params.gnomad_prepare},{params.maps_score} \
            scripts/maps_gcp.py  \
                -o $OUTPUT_GCP --intervals $(basename {input}) \
                --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
                --context_ht {config[gnomAD][context_ht]} \
                --mutation_ht {config[gnomAD][mutation_rate_ht]} \
                --genome_assembly {config[genome_assembly]} \
                --log $LOG_GCP

        echo download from gs://{config[bucket]}/${{OUTPUT_GCP}}
        gsutil cp -r "gs://{config[bucket]}/${{OUTPUT_GCP}}" {output}
        gsutil cp $LOG_GCP {log.local}
        """

rule prepare_gnomAD:
    threads: 10
    output:
        gnomAD_ht=directory(output_root/'gnomAD.ht')
    script: '../scripts/prepare_gnomad.py'

rule MAPS_local:
    input:
         intervals=rules.convert_bedfile.output,
         gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
          maps=directory(output_root/'MAPS/{variant_subset}_local.ht')
    log: str(output_root/'logs/hail_MAPS_{variant_subset}_local.log')
    threads: 10
    script: '../scripts/maps_score.py'
