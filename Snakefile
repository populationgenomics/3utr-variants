configfile: 'config.yml'
output_root = config['output_root']
config['assembly_ucsc'] = 'hg19' if config['genome_assembly'] == 'GRCh37' else 'hg38'

include: 'rules/download_data.smk'
include: 'rules/extract_PAS_features.smk'
include: 'rules/MAPS.smk'

rule all:
    input: rules.score_MAPS.output

rule dependency:
    output: output_root + '/dependency.svg'
    shell: 'snakemake --dag | dot -Tsvg -Grankdir=TB > {output}'
