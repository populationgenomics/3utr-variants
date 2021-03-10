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
    feature: pybedtools.Interval, signal: str, sequence: str, name: str = None
) -> Iterable[pybedtools.Interval]:
    """
    Create pybedtools.Interval objects of signal locations within given interval

    :param feature: interval in which to search for signal
    :param signal: sequence of signal motif
    :param sequence: sequence of feature interval
    :param name: value to overwrite in feature.name if specified
    :return: generator of signal locations
    """
    name = feature.name if name is None else name
    return (
        pybedtools.Interval(
            chrom=feature.chrom,
            start=feature.end - i.end()
            if feature.strand == '-'
            else feature.start + i.start(),
            end=feature.end - i.start()
            if feature.strand == '-'
            else feature.start + i.end(),
            name=name,
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

    :param feature: interval with interval.name containing annotation
        of form <PAS signal>|<conservation>
    :param sequence: sequence of feature interval
    :return: list of hexamer coordinates
    """
    signal_type, conservation = feature.name.split('|')  # pylint: disable=W0612

    if signal_type == 'NoPAS':
        return []
    if signal_type == 'OtherPAS':
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
    elif signal_type == 'Arich':
        signal_sequences = ['AAAAAA']
    else:
        signal_sequences = [signal_type]
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
    output = snakemake.output

    print('read database')
    df = pd.read_csv(pas_db, sep='\t')
    df['PSE'] = df['PSE'].str.rstrip('%').astype('float') / 100

    # create bedtools object/table
    df['Start'] = df['Position'] - 1
    df['Name'] = df['PAS Signal'] + '|' + df['Conservation']
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

    print('get flanking regions')
    bed_40nt = bed.slop(b=40, genome=genome)  # pylint: disable=unexpected-keyword-arg
    bed_100nt = bed.slop(b=100, genome=genome)  # pylint: disable=unexpected-keyword-arg

    print('extract hexamers')
    bed_40nt_us = bed.slop(  # pylint: disable=unexpected-keyword-arg
        l=40, r=0, s=True, genome=genome
    ).sequence(fi=fasta_file, s=True, fullHeader=True)
    sequences = FastaFile(bed_40nt_us.seqfn)

    hexamer_intervals = []
    for f in bed_40nt_us:
        seq = sequences.fetch(f'{f.chrom}:{f.start}-{f.stop}({f.strand})')
        hexamer_intervals.extend(get_hexamers(f, seq))
    hexamers_bed = pybedtools.BedTool(hexamer_intervals)

    # get stats
    with open(output.stats, 'w') as f:
        f.write(f'pA sites: {len(bed)}\n')
        f.write(f'hexamers: {len(hexamers_bed)}\n')

    print('save...')
    bed.saveas(output.PAS)
    bed_40nt.saveas(output.PAS_context_40nt)
    bed_100nt.saveas(output.PAS_context_100nt)
    hexamers_bed.saveas(output.PAS_hexamers)
