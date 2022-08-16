#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
	dx-download-all-inputs --parallel

    # Install packages from the python asset
	pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    echo "--------------Filtering and annotating excluded files -----------------"

    bedtools intersect -a $excluded_regions_path -b $panel_bed_path | sort | uniq > panel_exluded.bed
    head panel_exluded.bed
    # some panels may not be excluded so the panel_excluded.bed so have
    # an empty file outputted here. How to deal with the naming system?

    if [ -s panel_exluded.bed ]; then
        echo "Some panels are in excluded regions, annotation will be attempted."
        bedtools intersect -b $exons_hgnc_path -a panel_exluded.bed -wao > panel_excluded_genes.bed
        head panel_excluded_genes.bed
        python3 annotate_excluded_panel.py -e panel_excluded_genes.bed -p $panel_bed_path -r $excluded_regions_path -g $exons_gene_path
    else
		echo "Panel is not in excluded region"
        printf "Chrom\tStart\tEnd\tHGNC_ID\tTranscript\tExon\n" | tee touch $(echo ${excluded_regions_path##*/} | sed 's/.bed//g')_${panel_bed_path##*/}
	fi


    echo "--------------Outputting files -----------------"
    mkdir -p /home/dnanexus/out/annotated_excluded_file/

    out_file=$(find . -maxdepth 1 -name "00*.bed")
    mv $out_file /home/dnanexus/out/annotated_excluded_file/

    dx-upload-all-outputs
}
