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
