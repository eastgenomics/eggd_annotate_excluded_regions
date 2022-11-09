#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
	dx-download-all-inputs --parallel

    # Install packages from the python asset
	pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    # If no excluded regions were provided we can just create an empty annotate_regions file
    if [ "$excluded_regions" ] ; then
        echo "--------------Select first three columns of excluded region file -----------------"
        cut -f 1,2,3 $excluded_regions_path >  "$(basename $excluded_regions_path)"
        rm $excluded_regions_path
        mv $(basename $excluded_regions_path) in/excluded_regions/
        echo "--------------Filtering and annotating excluded files -----------------"
        # The cds exons file does not contain the extra regions, such as upstream
        # of exons and we need to annotate those if they are missing.
        # so we will have to concatonate the exons and additional regions files.

        # sometimes the additional file is not provided
        if [ "$additional_regions" ]; then
            # skip first row as its header than cat it to the end of the cds file
            sed 1d $additional_regions_path | cat $cds_hgnc_path - | sort -k1,1 -k2,2n > cds_exons.tsv
        else
            echo "No additional regions file is provided"
            mv $cds_hgnc_path cds_exons.tsv
        fi

        # sometimes a panel bed file is not provided
        if [ -z "$panel_bed" ]; then
            echo "No panel bed file provided, so the gCNV excluded regions is annotated."
            # -wao will keep regions that do and don't intersect
            bedtools intersect -b cds_exons.tsv -a $excluded_regions_path -wao > excluded_genes.bed
            python3 annotate_excluded_panel.py -e excluded_genes.bed -r $excluded_regions_path -c $cds_gene_path

        else
            echo "Panel bed file is provided, so the gCNV excluded regions will be intersected with panel bed file."
            # -wa will keep the a (excluded file) start and end rather than the start and end of both files
            bedtools intersect -a $excluded_regions_path -b $panel_bed_path -wa | sort | uniq > panel_excluded.bed
            # some panels may not be in the excluded file, so the panel_excluded.bed may be
            # empty. If its empty, the python script will error out. Therefore,
            # its easier to make an empty file with headers in the else statement. 
            if [ -s panel_excluded.bed ]; then
                echo "Some panel regions over lap with the gCNV excluded regions, these will be annotated."
                bedtools intersect -b cds_exons.tsv -a panel_excluded.bed -wao > panel_excluded_genes.bed
                python3 annotate_excluded_panel.py -e panel_excluded_genes.bed -p $panel_bed_path -r $excluded_regions_path -c $cds_gene_path
            else
                echo "Panel regions do not overlap with cnv calling excluded regions."
                printf "Chrom\tStart\tEnd\tLength\tGene_Symbol\tHGNC_ID\tTranscript\tExon\n" | tee touch annotated_$(echo ${excluded_regions_path##*/} | sed 's/.bed//g')_${panel_bed_path##*/}
            fi
        fi
    else
        echo "No excluded regions provided or excluded region empty"
        printf "Chrom\tStart\tEnd\tLength\tGene_Symbol\tHGNC_ID\tTranscript\tExon\n" | tee touch $(echo annotated_empty_excluded.bed)
    fi

    echo "--------------Outputting files -----------------"
    mkdir -p /home/dnanexus/out/annotated_excluded_file/

    out_file=$(ls | grep ^"annotated")
    echo $out_file
    mv $out_file /home/dnanexus/out/annotated_excluded_file/

    echo "--------------Converting bed to tsv -----------------"
    cd /home/dnanexus/out/annotated_excluded_file/
    linecount=$(wc -l $out_file  | cut -d " " -f 1)
    if [ "$linecount" -gt 1 ]; then
        echo "Adding 1bp";
        # check first row is header that starts with Chrom
        header=$(grep ^"Chrom" $out_file)
        if [ "$header" ]; then
            # get everything from the second line onwards as the first
            # line is the column names. Then add 1bp to start position. Then
            # join the first header to this file.
            tail -n +2 $out_file | sort -k1,1V -k2,2n | awk 'BEGIN {OFS="\t"}; {print $1,$2+1,$3,$4,$5,$6,$7,$8}' | cat <(head -n 1 $out_file ) - > ${out_file%.*}.tsv;
            # delete original annotated bed file
            rm $out_file;

        else
            echo "Incorrect output header"
        fi
    else
        echo "Empty file with headers only"
        mv $out_file ${out_file%.*}.tsv;
        rm $out_file;
    fi

    dx-upload-all-outputs
}
