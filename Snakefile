"""
Main Snakefile
Contains final rule, dependency graph rule, and connects Snakefiles from rules/
"""

from pathlib import Path

configfile: 'config.yml'
output_root = Path(config['output_root']).resolve()
config['assembly_ucsc'] = 'hg19' if config['genome_assembly'] == 'GRCh37' else 'hg38'

include: 'rules/download_data.smk'
include: 'rules/extract_UTR_features.smk'
include: 'rules/evaluation.smk'

rule all:
    input:
        expand(rules.gather_plot_single.input, run_location='local'),
        expand(
             rules.plots.output,
             run_location='local' if config['local'] else 'gcp'
         )

rule dependency:
    output:
        dag=multiext(str(output_root/'dependency_dag'), '.svg', '.png'),
        rulegraph=multiext(str(output_root/'dependency_rules'), '.svg', '.png')
    shell:
        """
        snakemake --dag | dot -Tsvg -Grankdir=TB > {output.dag[0]}
        snakemake --dag | dot -Tpng -Grankdir=TB > {output.dag[1]}
        snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {output.rulegraph[0]}
        snakemake --rulegraph | dot -Tpng -Grankdir=TB > {output.rulegraph[1]}
        """
