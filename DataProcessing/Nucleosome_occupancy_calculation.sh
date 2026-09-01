# Nucleosome occupancy calculation using DANPOS
# ATAC-seq and MNase-seq

cores=10

atac_bam_path="/path/to/ATAC/bam"
mnase_bam_path="/path/to/MNase/bam"
atac_danpos_path="/path/to/ATAC/danpos"
mnase_danpos_path="/path/to/MNase/danpos"
nucleosome_path="/path/to/nucleosome.bed"
genome_size="/path/to/genome_size.txt"

atac_samples="ATAC_sample_names"
mnase_samples="MNase_sample_names"


############################
# ATAC-seq
############################

for sample in ${atac_samples}; do

    mkdir -p ${atac_danpos_path}/${sample}/danpos

    # Select mononucleosomal fragments
    alignmentSieve \
        -b ${atac_bam_path}/${sample}_picard.bam \
        -p ${cores} \
        --ATACshift \
        --minFragmentLength 101 \
        --maxFragmentLength 215 \
        -o ${atac_danpos_path}/${sample}/${sample}_atac.bam

    # Sort and index BAM
    samtools sort -@ ${cores} \
        -o ${atac_danpos_path}/${sample}/${sample}_atac_sort.bam \
        ${atac_danpos_path}/${sample}/${sample}_atac.bam

    samtools index \
        ${atac_danpos_path}/${sample}/${sample}_atac_sort.bam

    # Call nucleosome positions and calculate occupancy using DANPOS
    python danpos.py dpos \
        ${atac_danpos_path}/${sample}/${sample}_atac_sort.bam \
        -o ${atac_danpos_path}/${sample}/danpos \
        -m 1 -q 1 -jd 73 -n N

    # Calculate normalization factor from total reads
    nf=$(echo "scale=10; 1000000/$(samtools view -@ ${cores} -c ${atac_bam_path}/${sample}_picard.bam)" | bc)

    # Normalize DANPOS occupancy signal to reads per million
    awk '/chrUn_/ {exit} {print}' \
        ${atac_danpos_path}/${sample}/danpos/pooled/${sample}_atac_sort.smooth.wig | \
    awk -v nf="$nf" '{
        if ($1 ~ /^-?[0-9]+(\.[0-9]+)?$/)
            printf "%.10f\n", $1*nf;
        else
            print $0
    }' \
    > ${atac_danpos_path}/${sample}/danpos/pooled/${sample}.wig

    # Convert normalized occupancy to bigWig
    wigToBigWig \
        ${atac_danpos_path}/${sample}/danpos/pooled/${sample}.wig \
        ${genome_size} \
        ${atac_danpos_path}/${sample}/danpos/pooled/${sample}.bw

done


############################
# MNase-seq
############################

for sample in ${mnase_samples}; do

    mkdir -p ${mnase_danpos_path}/${sample}

    # Call nucleosome positions and calculate normalized occupancy using DANPOS
    python danpos.py dpos \
        ${mnase_bam_path}/${sample}_sort.bam \
        -o ${mnase_danpos_path}/${sample} \
        -m 1 -q 1 -jd 73 -c 10000000

    # Convert normalized DANPOS occupancy signal to bigWig
    awk '/chrUn_/ {exit} {print}' \
        ${mnase_danpos_path}/${sample}/pooled/${sample}_sort.Fnor.smooth.wig \
        > ${mnase_danpos_path}/${sample}/pooled/${sample}.wig

    wigToBigWig \
        ${mnase_danpos_path}/${sample}/pooled/${sample}.wig \
        ${genome_size} \
        ${mnase_danpos_path}/${sample}/pooled/${sample}.bw

done
