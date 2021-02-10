bed_files = {
    'PolyADB_hexamers': rules.extract_PolyA_DB.output.PAS_hexamers,
    'Gencode_UTR': rules.extract_Gencode_UTR.output
}

rule convert_bedfile:
    input: lambda wildcards: bed_files[wildcards.variant_subset]
    output: output_root/'intervals/{variant_subset}_hail.txt'
    params:
        chr_style_hail = '' if config['genome_assembly'] == 'GRCh37' else 'chr'
    script: '../scripts/convert_bed.py'

rule MAPS_GCP:
    input: rules.convert_bedfile.output
    output: directory(output_root/'MAPS/{variant_subset}.ht')
    log: str(output_root/'logs/hail_MAPS_{variant_subset}.log')
    params:
        gnomad_prepare='../scripts/prepare_gnomAD.py',
        maps_score='../scripts/maps_score.py'
    shell:
        """
        OUTPUT=maps_{wildcards.variant_subset}.ht
        hailctl dataproc submit {config[cluster]} \
            --files {input} \
            --pyfiles {params.gnomad_prepare},{params.maps_score} \
            ../scripts/maps_gcp.py  -o $OUTPUT --intervals {input} \
            --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
            --context_ht {config[gnomAD][context_ht]} \
            --mutation_ht {config[gnomAD][mutation_rate_ht]} \
            --genome_assembly {config[genome_assembly]} \
            --subset_interval {config[gnomAD][subset]}

        #echo download from gs://{config[bucket]}/$OUTPUT
        #gsutil cp -r gs://{config[bucket]}/$OUTPUT {output}
        #gsutil cp gs://{config[bucket]}/logfile.log {log}
        """

rule prepare_gnomAD:
    threads: 10
    output:
        checkpoint=temp(directory(output_root/'gnomAD_tmp.ht')),
        gnomAD_ht=directory(output_root/'gnomAD.ht')
    script: '../scripts/prepare_gnomAD.py'

rule MAPS_local:
    input:
         intervals=rules.convert_bedfile.output,
         gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
          maps=directory(output_root/'MAPS/{variant_subset}_local.ht')
    log: str(output_root/'logs/hail_MAPS_{variant_subset}_local.log')
    threads: 10
    script: '../scripts/maps_score.py'
