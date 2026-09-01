# Cancer_type_enriched_nucleosome motif enrichment analysis

nucleosome_path=/path/to/cancer_type_enriched_nucleosome
genome=/path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
output_path=/path/to/homer/output
tmp_path=/path/to/homer/tmp
cancer_info=/path/to/cancer_info.txt

# Cancer type list
cancer_list=$(awk 'NR>1{print $4}' ${cancer_info} | sort -u | sed 's/TCGA-//')

for cancer in ${cancer_list}; do

    cancer_output=${output_path}/${cancer}

    # Create output directory
    mkdir -p ${cancer_output}

    # Extract cancer_type_enriched_nucleosome positions
    grep ${cancer} ${nucleosome_path} |
        awk '{print $1,$2,$3}' OFS="\t" \
        > ${tmp_path}/${cancer}.txt

    # Use all nucleosomes as the background set
    awk '{print $1,$2,$3}' OFS="\t" \
        ${nucleosome_path} \
        > ${tmp_path}/${cancer}_bg.txt

    # Identify enriched DNA sequence motifs using HOMER
    findMotifsGenome.pl \
        ${tmp_path}/${cancer}.txt \
        ${genome} \
        ${cancer_output} \
        -bg ${tmp_path}/${cancer}_bg.txt \
        -p 10 \
        -mset vertebrates

done
