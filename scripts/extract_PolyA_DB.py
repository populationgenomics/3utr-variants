import pandas as pd
import pybedtools

bed_file = str(snakemake.input[0])

df = pd.read_csv(bed_file, sep='\t')
df['PSE'] = df['PSE'].str.rstrip('%').astype('float') / 100
# filter dataset
# create bedtools object/table
df['Start'] = df['Position'] - 1
df['Name'] = df['Gene Symbol'] + '/' + df['Position'].map(str) + '/' + df['Strand']
bed_cols = [
    'Chromosome',  # chrom
    'Start',  # start
    'Position',  # end
    'Name',  # name
    'PSE',  # score
    'Strand'  # strand
]
bed = pybedtools.BedTool.from_dataframe(df[bed_cols])
# PAS only
bed.saveas(snakemake.output.PAS)
# PAS +/- 40nt
bed.slop(b=40, genome='hg19').saveas(snakemake.output.PAS_context_40nt)
# PAS +/- 100nt
bed.slop(b=100, genome='hg19').saveas(snakemake.output.PAS_context_100nt)
