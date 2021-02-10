# 3' UTR-variants

Identify variants occurring in 3’ UTRs that could disrupt the function of a transcript.

This project involves:

- retrieving reference genome using [genomepy](https://github.com/vanheeringen-lab/genomepy)
- extracing 3'UTR regions from gene annotation
- extracting PAS hexamers from [PolyA_DB3](https://exon.apps.wistar.org/polya_db/v3/misc/download.php)
- computing MAPS score on different 3'UTR variants locally and on the cloud using [hailctl dataproc](https://hail.is/docs/0.2/cloud/google_cloud.html)

## Pipeline

The workflow is implemented as a snakemake pipeline.
The dependencies are available as a conda environment.

```commandline
conda env create -f environment.yml
```

The `config.yml` file contains links to hail datasets and local files.
Local input files and output direcories can be modified accordingly.
In order to run the pipeline, simply call

```commandline
snakemake -n
```

for a dry run, or

```commandline
snakemake -j 10
```

to invoke the scripts.
