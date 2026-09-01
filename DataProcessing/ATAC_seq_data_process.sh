# ATAC-seq preprocessing and nucleosome positioning using DeNAPA and NucleoATAC

cores=10

fastq_path="/fastq/path/"
trim_path="/trimmed_fastq/path/"
bam_path="/bam/path/"
nucleosome_path="/nucleosome/path/"
genome_index="/path/to/GRCh38_index"
genome_fasta="/path/to/GRCh38.fa"
genome_size="/path/to/GRCh38.fa.fai"

# Sample list
list="ATAC_sample_names"

for sample in ${list}; do

    # Adapter and quality trimming
    trim_galore \
        --fastqc \
        -j ${cores} \
        --paired \
        -o ${trim_path} \
        ${fastq_path}/${sample}_1.fastq.gz \
        ${fastq_path}/${sample}_2.fastq.gz

    # Align paired-end reads to the GRCh38 reference genome
    bowtie2 \
        -x ${genome_index} \
        --very-sensitive \
        -X 2000 \
        -p ${cores} \
        -1 ${trim_path}/${sample}_1_val_1.fq.gz \
        -2 ${trim_path}/${sample}_2_val_2.fq.gz | \
    samtools sort \
        -@ ${cores} \
        -O BAM \
        -o ${bam_path}/${sample}_raw.bam -

    # Alignment quality assessment
    samtools flagstat \
        -@ ${cores} \
        ${bam_path}/${sample}_raw.bam \
        > ${bam_path}/qc/${sample}_raw.txt

    # Filter properly paired and high-quality mapped reads
    samtools view \
        -@ ${cores} \
        -bh \
        -f 2 \
        -q 30 \
        -F 2048 \
        ${bam_path}/${sample}_raw.bam \
        > ${bam_path}/${sample}_QC.bam

    # Remove PCR duplicates
    picard MarkDuplicates \
        INPUT=${bam_path}/${sample}_QC.bam \
        OUTPUT=${bam_path}/${sample}_picard.bam \
        METRICS_FILE=${bam_path}/qc/${sample}.metrics \
        REMOVE_DUPLICATES=true

    # Sort and index BAM
    samtools sort \
        -@ ${cores} \
        -o ${bam_path}/${sample}_sort.bam \
        ${bam_path}/${sample}_picard.bam

    samtools index \
        -@ ${cores} \
        ${bam_path}/${sample}_sort.bam


    # Call nucleosome positions using DeNAPA
    mkdir -p ${nucleosome_path}/denopa/${sample}

    denopa \
        -i ${bam_path}/${sample}_picard.bam \
        -o ${nucleosome_path}/denopa/${sample} \
        --proc ${cores} \
        -c "chr[1-9][0-9]{,1}|chrX" \
        -r


    # Generate broad accessible regions for NucleoATAC
    mkdir -p ${nucleosome_path}/nucleoatac/${sample}/Macs2_out

    macs2 callpeak \
        -t ${bam_path}/${sample}_picard.bam \
        -n ${sample} \
        --shift -75 \
        --extsize 150 \
        --nomodel \
        --nolambda \
        --keep-dup all \
        -g 2747877702 \
        --broad \
        --broad-cutoff 0.01 \
        --outdir ${nucleosome_path}/nucleoatac/${sample}/Macs2_out

    # Extend and merge candidate regions
    bedtools slop \
        -i ${nucleosome_path}/nucleoatac/${sample}/Macs2_out/${sample}_peaks.broadPeak \
        -g ${genome_size} \
        -b 1000 | \
    bedtools merge -i - \
        > ${nucleosome_path}/nucleoatac/${sample}/Macs2_out/${sample}_slop_peaks.broadPeak


    # Call nucleosome positions using NucleoATAC
    cd ${nucleosome_path}/nucleoatac/${sample}

    nucleoatac run \
        --bed ${nucleosome_path}/nucleoatac/${sample}/Macs2_out/${sample}_slop_peaks.broadPeak \
        --bam ${bam_path}/${sample}_picard.bam \
        --fasta ${genome_fasta} \
        --out ${sample}_slop \
        --cores ${cores}

done
