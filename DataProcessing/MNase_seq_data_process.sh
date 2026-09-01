# MNase-seq preprocessing and nucleosome positioning analysis using DANPOS

cores=10

qc_path=/qc/path/
fastq_path=/fastq/path/
trim_path=/trimmed_fastq/path/
bam_path=/bam/path/
danpos_path=/danpos/path/

# Sample list
list="MNase-seq_sample_names"

for sample in ${list}; do

    # Adapter and quality trimming
    trim_galore --fastqc -j ${cores} \
        --paired \
        -o ${trim_path} \
        ${fastq_path}/${sample}_1.fastq.gz \
        ${fastq_path}/${sample}_2.fastq.gz

    # Align reads to the human reference genome using Bowtie2
    bowtie2 -x GRCh38_index \
        --very-sensitive -X 2000 -p ${cores} \
        -1 ${trim_path}/${sample}_1_val_1.fq.gz \
        -2 ${trim_path}/${sample}_2_val_2.fq.gz | \
        samtools sort -@ ${cores} -O bam \
        -o ${bam_path}/${sample}_raw.bam -

    # Retain properly paired, high-quality alignments
    samtools view -@ ${cores} -bh \
        -f 2 -q 30 -F 2048 \
        ${bam_path}/${sample}_raw.bam \
        > ${bam_path}/${sample}_QC.bam


    # Remove PCR duplicates
    picard MarkDuplicates \
        INPUT=${bam_path}/${sample}_QC.bam \
        OUTPUT=${bam_path}/${sample}_picard.bam \
        METRICS_FILE=${qc_path}/${sample}.mat \
        REMOVE_DUPLICATES=true

    # Sort and index BAM files
    samtools sort -@ ${cores} \
        -o ${bam_path}/${sample}_sort.bam \
        ${bam_path}/${sample}_picard.bam

    samtools index -@ ${cores} \
        ${bam_path}/${sample}_sort.bam

    # Nucleosome positioning analysis using DANPOS
    mkdir -p ${danpos_path}/${sample}

    python danpos.py dpos \
        ${bam_path}/${sample}_sort.bam \
        -o ${danpos_path}/${sample} \
        -m 1 -n F

done
