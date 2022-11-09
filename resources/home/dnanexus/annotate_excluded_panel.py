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
        help='Excluded regions file intersected with requested panel regions',
        required=True
        )

    parser.add_argument(
        '-p', '--panel',
        help='Panel bed file',
        required=False
        )

    parser.add_argument(
        '-r', '--excluded_region',
        help='Excluded regions file',
        required=True
        )

    parser.add_argument(
        '-c', '--cds',
        help='cds file by gene symbols',
        required=True
        )

    args = parser.parse_args()

    return args


def read_data(args):
    """Reads in data from the arguement inputs. Also checks whether panel
    bed file is provided.

    Args:
        args (object): parse_args object containing all arg input info

    Returns:
        exc_panel (pd datatframe): gCNV exluded file (either interesected
                                    eith panel or not)
        panel (pd dataframe OR null): panel file
        cds_gene (pd dataframe): cds file containing gene symbols
    """

    # read data in
    exc_panel = pd.read_csv(args.excluded_panel, sep="\t", header=None)
    cds_gene = pd.read_csv(args.cds, sep="\t",
                            header=None, dtype='unicode')

    # Check input files have expected columns and read in data
    exc_panel_col_names = ["chr_exluded", "pos_start_excluded",
                            "pos_end_excluded",
                            "chr_GCF", "pos_start_GCF", "pos_end_GCF",
                            "HGNC_ID", "transcript", "exon", "num"]
    if len(exc_panel.columns) != len(exc_panel_col_names):
        raise Exception("excluded_panel file '{}' does not "
                        "contain expected columns".format(
                            args.excluded_panel
                            ))
    exc_panel.columns = exc_panel_col_names

    cds_gene_col_names = ["Chr", "Start",
                        "End", "Gene_Symbol", "Transcript",
                        "Exon"]
    if len(cds_gene.columns) != len(cds_gene_col_names):
        raise Exception("cds_gene file does not "
                        "contain expected columns")
    cds_gene.columns = cds_gene_col_names

    if args.panel is None:
        print('No panel bed inputted into python script.')
        panel  = ""
    else:
        panel = pd.read_csv(args.panel, sep="\t", header=None)
        # Check input files have expected columns and read in data
        panel_col_names = ["chr", "pos_start", "pos_end", "transcript"]
        if len(panel.columns) != len(panel_col_names):
            raise Exception(
                "Panel file '{}' does not contain "
                "expected columns".format(args.panel))
        panel.columns = panel_col_names


    return exc_panel, panel, cds_gene


def main():

    args = parse_args()

    exc_panel, panel, cds_gene = read_data(args)

    if args.panel is not None:
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
    else:
        exc_panel_transcript = exc_panel

    # select excluded columns, HGNCID, transcript, exon
    exc_panel_transcript_subset = exc_panel_transcript[[
                                    "chr_exluded", "pos_start_excluded",
                                    "pos_end_excluded", "HGNC_ID",
                                    "transcript", "exon"]]
    # rename columns
    exc_panel_transcript_subset.columns = ["Chrom", "Start", "End",
                                        "HGNC_ID", "Transcript", "Exon"]
    # Calculate length of annotated excluded region
    length = exc_panel_transcript_subset.loc[:, 'End'] - exc_panel_transcript_subset.loc[:, 'Start']
    exc_panel_transcript_subset.insert(6, "Length", length, True)

    # lets add the gene symbol now
    # take the gene & transcript info
    cds_gene_subset = cds_gene[["Gene_Symbol", "Transcript"]]
    cds_gene_subset = cds_gene_subset.drop_duplicates()

    # left join on transcript to get the gene symbol
    df = exc_panel_transcript_subset.merge(cds_gene_subset,
                                            on='Transcript', how='left')
    df2 = df.replace(np.nan, '.', regex=True)

    # reorder columns
    df2 = df2[["Chrom", "Start", "End", "Length",
                "Gene_Symbol", "HGNC_ID","Transcript", "Exon"]]

    # name output file excluded_file + panel
    excluded_name = args.excluded_region
    excluded_name = excluded_name.split("/")[-1].rstrip(".bed")

    if args.panel is not None:
        panel_name = args.panel
        panel_name = panel_name.split("/")[-1].rstrip(".bed")
        output_filename = "annotated_" + excluded_name + "_" + panel_name + ".bed"
    else:
        output_filename = "annotated_" + excluded_name  + ".bed"

    df2.to_csv(
            output_filename,
            sep="\t", index=False, header=True
            )


if __name__ == "__main__":

    main()
