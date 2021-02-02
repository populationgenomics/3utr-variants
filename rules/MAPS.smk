output_root = config['output_root']

rule prepare_gnomAD:
    threads: 10
    output:
        checkpoint=temp(directory(output_root + '/gnomAD_tmp.ht')),
        gnomAD_ht=directory(output_root + '/gnomAD.ht')
    script: '../scripts/prepare_gnomAD.py'

rule score_MAPS:
    # TODO: generalise by DB type
    input:
         bed=rules.extract_Gencode_UTR.output,
         gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
          maps=directory(output_root + '/MAPS/Gencode.ht')
    log: output_root + '/logs/hail_MAPS_Gencode.log'
    threads: 10
    script: '../scripts/maps_score.py'
