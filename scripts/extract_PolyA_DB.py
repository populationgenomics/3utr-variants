import pandas as pd
import pybedtools
from pysam import FastaFile
import re


def create_signal_interval(feature, signal, sequence):
    return (
        pybedtools.Interval(
            chrom=feature.chrom,
            start=feature.start + i.start(),
            end=feature.start + i.start() + 6,
            name=signal,
            score=feature.score,
            strand=feature.strand
        )
        for i in re.finditer(signal, sequence, re.IGNORECASE)
    )


def get_hexamers(feature: pybedtools.Interval, sequence: str):
    if feature.name == "NoPAS":
        return []
    elif feature.name == "Other PAS":
        signal_sequences = ['AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA', 'AATATA', 'AATACA', 'AATAGA', 'AAAAAG', 'ACTAAA']
    elif feature.name == "Arich":
        signal_sequences = ['AAAAAA']
    else:
        signal_sequences = [feature.name]
    dna_signals = [s.replace('U', 'T') for s in signal_sequences]

    hexamer_list = []
    for signal in dna_signals:
        intervals = create_signal_interval(feature, signal, sequence)
        hexamer_list.extend(intervals)
    return hexamer_list


if __name__ == '__main__':
    genome = snakemake.config['assembly_ucsc']
    pas_db = str(snakemake.input.db)
    fasta_file = str(snakemake.input.fasta)

    out_pas = snakemake.output.PAS
    out_40nt = snakemake.output.PAS_context_40nt
    out_100nt = snakemake.output.PAS_context_100nt
    out_hex = snakemake.output.PAS_hexamers

    print('read database')
    df = pd.read_csv(pas_db, sep='\t')
    df['PSE'] = df['PSE'].str.rstrip('%').astype('float') / 100
    # TODO: filter dataset

    # create bedtools object/table
    df['Start'] = df['Position'] - 1
    df['Name'] = df['PAS Signal']  # df['Gene Symbol']  + '/' + df['Position'].map(str) + '/' + df['Strand']
    bed_cols = [
        'Chromosome',  # chrom
        'Start',  # start
        'Position',  # end
        'Name',  # name
        'PSE',  # score
        'Strand'  # strand
    ]
    bed = pybedtools.BedTool.from_dataframe(df[bed_cols].drop_duplicates())
    print('BedTools object created')

    # extract hexamers form 40nt upstream
    bed_40nt = bed.slop(l=40, r=0, s=True, genome=genome).sequence(fi=fasta_file, s=True, fullHeader=True)
    sequences = FastaFile(bed_40nt.seqfn)

    print('extract hexamers')
    hexamer_intervals = []
    for feature in bed_40nt:
        sequence = sequences.fetch(f"{feature.chrom}:{feature.start}-{feature.stop}({feature.strand})")
        hexamer_intervals.extend(get_hexamers(feature, sequence))

    print('save...')
    bed.saveas(out_pas)
    bed.slop(b=40, genome=genome).saveas(out_40nt)
    bed.slop(b=100, genome=genome).saveas(out_100nt)
    pybedtools.BedTool(hexamer_intervals).sort().merge().saveas(out_hex)
