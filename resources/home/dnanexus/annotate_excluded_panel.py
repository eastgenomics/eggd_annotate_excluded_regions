#!/usr/bin/python


import argparse
import numpy as np
import pandas as pd

def parse_args():
    """Parse through arguments
    Returns:
        args: Variable that you can extract relevant
        arguments inputs needed
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

    parser.add_argument(
        '-r', '--excluded_region',
        help='Excluded region file',
        required=True
        )

    parser.add_argument(
        '-g', '--exons',
        help='GCF exons file by gene symbols',
        required=True
        )

    args = parser.parse_args()

    return args


def main():

    args = parse_args()

    # read data in
    exc_panel = pd.read_csv(args.excluded_panel, sep="\t", header=None)
    panel = pd.read_csv(args.panel, sep="\t", header=None)
    exons_gene = pd.read_csv(args.exons, sep="\t", header=None)

    # Check input files have expected columns and read in data
    exc_panel_col_names = ["chr_exluded", "pos_start_excluded",
                            "pos_end_excluded", "strand", "dot",
                            "chr_GCF", "pos_start_GCF", "pos_end_GCF",
                            "HGNC_ID", "transcript", "exon", "num"]
    if len(exc_panel.columns) != len(exc_panel_col_names):
        raise Exception("excluded_panel file '{}' does not "
                        "contain expected columns".format(
                            args.excluded_panel
                            ))
    exc_panel.columns = exc_panel_col_names

    panel_col_names = ["chr", "pos_start", "pos_end", "transcript"]
    if len(panel.columns) != len(panel_col_names):
        raise Exception(
            "Panel file '{}' does not contain "
            "expected columns".format(args.panel))
    panel.columns = panel_col_names

    # Check both files have expected columns and read in data
    exons_gene_col_names = ["Chr", "Start",
                        "End", "Gene_Symbol", "Transcript",
                        "Exon"]
    if len(exons_gene.columns) != len(exons_gene_col_names):
        raise Exception("exons_gene file does not "
                        "contain expected columns")
    exons_gene.columns = exons_gene_col_names

    # get the transcripts in panel
    panel_transcripts = list(panel['transcript'].unique())
    # Add the dot for cases where its not exonic so they have a
    # dot instead of a transcript.
    panel_transcripts.append(".")

    # keep rows that have panel transcript in the exc_panel as exc_panel
    # has many transcript to gene
    exc_panel_transcript = exc_panel.loc[
                        exc_panel["transcript"].isin(panel_transcripts)
                        ]

    # select excluded columns, HGNCID, transcript, exon
    exc_panel_transcript_subset = exc_panel_transcript[[
                                    "chr_exluded", "pos_start_excluded",
                                    "pos_end_excluded", "HGNC_ID",
                                    "transcript", "exon"]]
    # rename columns
    exc_panel_transcript_subset.columns = ["Chrom", "Start", "End",
                                        "HGNC_ID", "Transcript", "Exon"]
    # Calculate length of annotated excluded region
    exc_panel_transcript_subset['Length'] = exc_panel_transcript_subset['End'] - exc_panel_transcript_subset['Start']

    # lets add the gene symbol now
    # take the gene & transcript info
    exons_gene_subset = exons_gene[["Gene_Symbol", "Transcript"]]
    exons_gene_subset = exons_gene_subset.drop_duplicates()

    # left join on transcript to get the gene symbol
    df = exc_panel_transcript_subset.merge(exons_gene_subset,
                                            on='Transcript', how='left')
    df2 = df.replace(np.nan, '.', regex=True)

    # reorder columns
    df2 = df2[["Chrom", "Start", "End", "Length",
                "Gene_Symbol", "HGNC_ID","Transcript", "Exon"]]

    # name output file excluded_file + panel
    panel_name = args.panel
    panel_name = panel_name.split("/")[-1]
    panel_name = panel_name.split(".bed")[0]

    excluded_name = args.excluded_region
    excluded_name = excluded_name.split("/")[-1]
    excluded_name = excluded_name.split(".bed")[0]

    output_filename = excluded_name + "_" + panel_name + ".bed"
    print(output_filename)

    df2.to_csv(
            output_filename,
            sep="\t", index=False, header=True
            )


if __name__ == "__main__":

    main()
