"""
Extract PolyA sites and cis-elements & flanking regions from PolyA_DB3

When run as script:
    Input:
        PolyA_DB database file
        reference sequence in FASTA format
    Output: bed files (see rules/extract_UTR_features.smk)
"""

import re
from typing import Iterable
from typing import List

import pandas as pd
import pybedtools
from pysam import FastaFile  # pylint: disable=no-name-in-module


def create_signal_interval(
    feature: pybedtools.Interval, signal: str, sequence: str
) -> Iterable[pybedtools.Interval]:
    """
    Create pybedtools.Interval objects of signal locations within given interval

    :param feature: interval in which to search for signal
    :param signal: sequence of signal motif
    :param sequence: sequence of feature interval
    :return: generator of signal locations
    """
    return (
        pybedtools.Interval(
            chrom=feature.chrom,
            start=feature.end - i.end()
            if feature.strand == '-'
            else feature.start + i.start(),
            end=feature.end - i.start()
            if feature.strand == '-'
            else feature.start + i.end(),
            name=signal,
            score=feature.score,
            strand=feature.strand,
        )
        for i in re.finditer(signal, sequence, re.IGNORECASE)
    )


def get_hexamers(
    feature: pybedtools.Interval, sequence: str
) -> List[pybedtools.Interval]:
    """
    Retrieve all hexamers of regions from within given interval

    :param feature: interval of interest
    :param sequence: sequence of feature interval
    :return: list of hexamer coordinates
    """
    if feature.name == 'NoPAS':
        return []
    if feature.name == 'OtherPAS':
        signal_sequences = [
            'AGTAAA',
            'TATAAA',
            'CATAAA',
            'GATAAA',
            'AATATA',
            'AATACA',
            'AATAGA',
            'AAAAAG',
            'ACTAAA',
        ]
    elif feature.name == 'Arich':
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
    pas_db = snakemake.input.db.__str__()
    fasta_file = snakemake.input.fasta.__str__()

    print('read database')
    df = pd.read_csv(pas_db, sep='\t')
    df['PSE'] = df['PSE'].str.rstrip('%').astype('float') / 100
    if snakemake.wildcards['filter'] == 'conserved':
        df = df[df['Conservation'] == 'Yes']

    # create bedtools object/table
    df['Start'] = df['Position'] - 1
    df['Name'] = df[
        'PAS Signal'
    ]  # df['Gene Symbol']  + '/' + df['Position'].map(str) + '/' + df['Strand']
    bed_cols = [
        'Chromosome',  # chrom
        'Start',  # start
        'Position',  # end
        'Name',  # name
        'PSE',  # score
        'Strand',  # strand
    ]
    bed = pybedtools.BedTool.from_dataframe(df[bed_cols].drop_duplicates())
    print('BedTools object created')

    print('extract hexamers')
    bed_40nt = bed.slop(  # pylint: disable=unexpected-keyword-arg
        l=40, r=0, s=True, genome=genome
    ).sequence(fi=fasta_file, s=True, fullHeader=True)
    sequences = FastaFile(bed_40nt.seqfn)

    hexamer_intervals = []
    for f in bed_40nt:
        seq = sequences.fetch(f'{f.chrom}:{f.start}-{f.stop}({f.strand})')
        hexamer_intervals.extend(get_hexamers(f, seq))

    print('save...')
    out = snakemake.output
    # fmt: off
    bed.saveas(out.PAS)
    # pylint: disable=unexpected-keyword-arg
    bed.slop(b=40, genome=genome).saveas(out.PAS_context_40nt)
    # pylint: disable=unexpected-keyword-arg
    bed.slop(b=100, genome=genome).saveas(out.PAS_context_100nt)
    pybedtools.BedTool(hexamer_intervals).saveas(out.PAS_hexamers)

    with open(out.stats, 'w') as f:
        f.write(f'pA sites: {len(bed)}\n')
        f.write(f'hexamers: {len(hexamer_intervals)}\n')
