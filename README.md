# eggd_annotate_excluded_regions

## What does this app do?

This app takes in the excluded region file generated by gCNV_call for every run to show regions that CNV calling wasn't performed in. The excluded regions file will be filtered for the panel, then annotated for the HGNC ID, gene symbol and exon number using the exons file. The length of the annotated excluded regions will be calculated as well.

## What are the inputs?
- Excluded regions file
- Panel bed
- Exons HGNC ID file
- Exons gene symbols file
- Additional regions file (chromosome, start, end, gene_panel, transcript, exon)

## What are the outputs?
- Annotated excluded region filtered for panel regions.

## Where is this app applicable?
For germline CNV reports.

### This app was made by EMEE GLH