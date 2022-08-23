#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
	dx-download-all-inputs --parallel

    # Install packages from the python asset
	pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    echo "--------------Filtering and annotating excluded files -----------------"
    # sometimes a panel bed file is not provided
    if [ -z "$panel_bed" ]; then
        echo "No panel bed file provided, so the gCNV excluded regions is annotated."
        bedtools intersect -b $cds_hgnc_path -a $excluded_regions_path -wao > excluded_genes.bed
        head excluded_genes.bed
        python3 annotate_excluded_panel.py -e excluded_genes.bed -r $excluded_regions_path -c $cds_gene_path

    else
        echo "Panel bed file is provided, so the gCNV excluded regions will be interested with panel bed file."
        bedtools intersect -a $excluded_regions_path -b $panel_bed_path | sort | uniq > panel_excluded.bed
        head panel_excluded.bed
        # some panels may not be excluded so the panel_excluded.bed so have
        # an empty file outputted here.
        if [ -s panel_excluded.bed ]; then
            echo "Some panel regions over lap with the gCNV excluded regions, these will be annotated."
            bedtools intersect -b $cds_hgnc_path -a panel_excluded.bed -wao > panel_excluded_genes.bed
            head panel_excluded_genes.bed
            python3 annotate_excluded_panel.py -e panel_excluded_genes.bed -p $panel_bed_path -r $excluded_regions_path -c $cds_gene_path
        else
            echo "Panel regions do not overlap with cnv calling excluded regions."
            printf "Chrom\tStart\tEnd\tHGNC_ID\tTranscript\tExon\n" | tee touch $(echo ${excluded_regions_path##*/} | sed 's/.bed//g')_${panel_bed_path##*/}
        fi
    fi

    echo "--------------Outputting files -----------------"
    mkdir -p /home/dnanexus/out/annotated_excluded_file/

    out_file=$(find . -maxdepth 1 -name "00*.bed")
    mv $out_file /home/dnanexus/out/annotated_excluded_file/

    dx-upload-all-outputs
}
