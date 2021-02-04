from pathlib import Path

configfile: 'config.yml'
output_root = Path(config['output_root']).resolve()
config['assembly_ucsc'] = 'hg19' if config['genome_assembly'] == 'GRCh37' else 'hg38'

include: 'rules/download_data.smk'
include: 'rules/extract_PAS_features.smk'
include: 'rules/evaluation.smk'

rule all:
    input:
        expand(
            rules.MAPS_local.output if config['local'] else rules.MAPS_GCP.output,
            variant_subset=[
                'PolyADB_hexamers',
                #'Gencode_UTR'
            ]
        )

rule dependency:
    output: output_root/'dependency.svg'
    shell: 'snakemake --dag | dot -Tsvg -Grankdir=TB > {output}'
