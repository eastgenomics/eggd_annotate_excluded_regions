#!/usr/bin/python

import pandas as pd
import argparse


def parse_args():
    """Parse through arguements
    Returns:
        args: Variable that you can extract relevant
        arguements inputs needed
    """
    # Read in arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-e', '--excluded_panel',
        help='Excluded file filtered to panel',
        required=True
        )

    parser.add_argument(
        '-p', '--panel',
        help='Panel bed file',
        required=True
        )

    args = parser.parse_args()

    return args



def main():

    args = parse_args()

    # read data in
    exc_panel = pd.read_csv(args.excluded_panel, sep="\t", header=None)
    panel  = pd.read_csv(args.panel, sep="\t", header=None)

    # Check both files have expected columns and read in data
    exc_panel_col_names = ["chr_exluded", "pos_start_excluded", "pos_end_excluded", "strand", "dot", "chr_GCF", "pos_start_GCF", "pos_end_GCF", "HGNC_ID", "transcript", "exon"]
    if len(exc_panel.columns) != len(exc_panel_col_names):
        raise Exception("excluded_panel file '{}' does not "
        "contain expected columns".format(args.excluded_panel))
    exc_panel.columns = exc_panel_col_names

    panel_col_names = ["chr", "pos_start", "pos_end", "transcript"]
    if len(panel.columns) != len(panel_col_names):
        raise Exception("Panel file '{}' does not contain "
        "expected columns".format(args.panel))
    panel.columns = panel_col_names

    # get the transcripts in panel
    panel_transcripts = list(panel['transcript'].unique())

    # keep rows that have panel transcript in the exc_panel as exc_panel
    # has many transcript to gene
    exc_panel_transcript = exc_panel.loc[exc_panel["transcript"].isin(panel_transcripts) ]

    # select excluded columns, HGNCID, transcript, exon
    exc_panel_transcript_subset = exc_panel_transcript[["chr_exluded", "pos_start_excluded","pos_end_excluded",  "HGNC_ID", "transcript", "exon"]]
    # rename columns
    exc_panel_transcript_subset.columns = ["Chr", "Pos_start","Pos_end",  "HGNC_ID", "Transcript", "Exon"]

    exc_panel_transcript_subset.to_csv(
            'annotated_excluded_panel_region.bed',
            sep="\t", index=False, header=True
            )

if __name__ == "__main__":

    main()

