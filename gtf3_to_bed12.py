# Convert a gtf file to a bed12+ file
# Assumes standard gtf as per https://useast.ensembl.org/info/website/upload/gff.html
# Output is standard bed12 format as per https://bedtools.readthedocs.io/en/latest/content/general-usage.html
# Usage: python gtf3_to_bed12.py <gtf_file> <output_file>
from tqdm import tqdm
from p_tqdm import p_map
import pandas as pd
import argparse
import sys
import numpy as np


def transcript_rows_to_bed_line(g, bed12plus=False):
    """
    Function to convert a group of rows from a gtf file to a bed12+ line. Used with df.groupby().apply()

    params 
    g: grouped pandas DataFrame
        Must be a gtf file with named columns using standard gtf format. DataFrame must have been grouped with df.groupby()
    bed12plus: bool
        If true, keep transcipt_id as 13th column. Default is False. Current not implemented.

    returns: pandas Series
        A single row of bed12+ format

    example:
    gtf = pd.read_csv('genode_file.gtf', sep='\t', comment='#', header=gtf_cols)
    bed = gtf.groupby('transcript_id').apply(transcript_rows_to_bed_line)
    """
    if bed12plus:
        raise NotImplementedError("Sorry, bed12plus is not yet implemented")

    # get general block information where feature is transcript and blocks are exons
    start, end = g.loc[g['feature'] == 'transcript',
                       ['start', 'end']].to_numpy()[0]
    blocks = g.loc[g['feature'] == 'exon', ['start', 'end']].assign(
        size=lambda x: x['end'] - x['start'])
    strand = g['strand'].iat[0].strip()

    def reverse_if_negative_strand(s):
        if strand == '+':
            return s
        else:
            return s[::-1]

    bed_line = pd.Series({
        'chrom':
        g['seqname'].iat[0],
        'start':
        start - 1,
        'end':
        end - 1,
        'name':
        g['gene_name'].iat[0],
        'score':
        0,
        'strand':
        strand,
        'thickstart':
        start - 1,
        'thickend':
        end - 1,
        'itemRGB':
        0,
        'blockCount':
        blocks.shape[0],
        'blocksizes':
        ','.join(reverse_if_negative_strand(blocks['size'].astype(str))),
        'blockstarts':
        ','.join(
            reverse_if_negative_strand(
                ((blocks['start']) - start).astype(str)))
    })

    if bed12plus:
        bed_line['transcript_id'] = g['transcript_id'].iat[0]

    return bed_line


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Convert gtf file to bed12+ format')
    argparser.add_argument('-g', '--gtf_file', type=str, help='Input gtf file')
    argparser.add_argument(
        '-b',
        '--bed12plus',
        default='False',
        help=
        'If true, keep transcript ID as additional column. If False, output is standrd bed12 format. Not implemented yet'
    )
    args = argparser.parse_args()

    tqdm.pandas(desc="Processing GTF transcripts")
    in_file = args.gtf_file
    bed12plus = True if args.bed12plus.lower() == 'true' else False

    gtf_cols = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
        'frame', 'attribute'
    ]
    dtypes = dict(
        zip(gtf_cols,
            ['str', 'str', 'str', 'int', 'int', 'str', 'str', 'str', 'str']))
    df = pd.read_csv(in_file,
                     sep='\t',
                     comment='#',
                     names=gtf_cols,
                     dtype=dtypes)

    # by default, we chunk the gtf file into 1,000 chunks
    chunks = np.array_split(df, 1000)
    processed = []

    def parse_attribute(chunk):
        return chunk.assign(transcript_id=chunk.iloc[:, 8].str.extract(
            '.*transcript_id "(ENSM?U?S?T[0-9\.]+)";'
        )).assign(gene_id=chunk.iloc[:, 8].str.extract(
            '.*gene_id "(ENSM?U?S?G[0-9\.]+)";')).assign(
                gene_name=chunk.iloc[:, 8].str.extract('.*gene_name "(\w+)";'))
    
    # This part here parses the attribute field, and is parallelized
    processed = p_map(parse_attribute, chunks,desc='Parsing of attributes column, 1,000 chunks')
    df = pd.concat(processed)

    df['gene_name'] = df['gene_name'].fillna(df['gene_id']).fillna(
        df['transcript_id'])
#    pd.concat(p_map(lambda x: transcript_rows_to_bed_line(x[1], bed12plus=bed12plus),
#              df.groupby('transcript_id')), axis=1).T.to_csv(sys.stdout, sep='\t', index=False, header=False)
    #  Here progress_apply is not parallelized, but parallelizing over a groupby doesn't seem to be efficient in practice (as implement in the line above)
    df.groupby('transcript_id').progress_apply(
        transcript_rows_to_bed_line, bed12plus=bed12plus).to_csv(sys.stdout, sep='\t', index=False, header=False)
