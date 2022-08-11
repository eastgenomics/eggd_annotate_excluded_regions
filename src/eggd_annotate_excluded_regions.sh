#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
	dx-download-all-inputs --parallel

    # Install packages from the python asset
	pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    echo "--------------Filtering and annotating excluded files -----------------"

    bedtools intersect -a $excluded_regions_path -b $panel_bed_path > panel_exluded.bed

    bedtools intersect -b $exons_file_path -a panel_exluded.bed -wa -wb > panel_excluded_genes.bed

    python3 annotate_excluded_panel.py -e panel_excluded_genes.bed -p $panel_bed_path -r $excluded_regions_path

    echo "--------------Outputting files -----------------"
    mkdir -p /home/dnanexus/out/annotated_excluded_file/

    out_file=$(find . -maxdepth 1 -name "00*.bed")
    mv $out_file /home/dnanexus/out/annotated_excluded_file/

    dx-upload-all-outputs
}
