#!/bin/bash
set -eo pipefail
SECONDS=0
cd ~/RNAseq_pipeline/ || { echo "Failed to enter pipeline directory"; exit 1; }

# Configuration
THREADS=10
STRANDEDNESS_FILE="strandedness.txt"
RAW_DIR="data/raw_fastq"
TRIMMED_DIR="data/trimmed_fastq"
QC_DIR="data/QC_trimmed_fastq"
ALIGN_DIR="HISAT2"
QUANT_DIR="quants"
SALMON_INDEX="salmon_index"
HISAT2_INDEX="HISAT2/grcm38/genome"
ANNOTATION="annotation.gtf"

# Create directories
mkdir -p "$TRIMMED_DIR" "$QC_DIR" "$ALIGN_DIR" "$QUANT_DIR" || { 
    echo "Directory creation failed"; exit 1; 
}

# Validation functions
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Missing required file $1" >&2
        exit 1
    fi
}

check_dir_exists() {
    if [ ! -d "$1" ]; then
        echo "ERROR: Directory $1 does not exist" >&2
        exit 1
    fi
}

check_dir_exists "$RAW_DIR"
check_file_exists "$ANNOTATION"
check_file_exists "$STRANDEDNESS_FILE"

# Strandedness function
get_strandedness() {
    local sample=$1
    local tool=$2

    strandedness=$(awk -v samp="$sample" '$1 == samp {print $2}' "$STRANDEDNESS_FILE")
    if [ -z "$strandedness" ]; then
        echo "ERROR: Sample $sample not found in $STRANDEDNESS_FILE" >&2
        exit 1
    fi

    case $tool in
        "HISAT2")
            case "$strandedness" in
                "SR") echo "R" ;;
                "SF") echo "F" ;;
                *)    echo "" ;;
            esac ;;
        "featureCounts")
            case "$strandedness" in
                "SR") echo 2 ;;
                "SF") echo 1 ;;
                *)    echo 0 ;;
            esac ;;
        "HTSeq")
            case "$strandedness" in
                "SR") echo "reverse" ;;
                "SF") echo "yes" ;;
                *)    echo "no" ;;
            esac ;;
        *)
            echo "ERROR: Unknown tool $tool" >&2
            exit 1 ;;
    esac
}

# Validate FASTQ files with _R1_ pattern
if [ -z "$(ls -A "$RAW_DIR"/*_R1_*.fastq.gz 2>/dev/null)" ]; then
    echo "ERROR: No R1 FASTQ files found in $RAW_DIR" >&2
    exit 1
fi

# STEP 1: Paired-end trimming with corrected filename handling
for r1_file in "$RAW_DIR"/*_R1_*.fastq.gz; do
    # Derive R2 filename from R1 filename
    r2_file="${r1_file/_R1_/_R2_}"
    check_file_exists "$r2_file"
    
    # Extract base name without R1 suffix
    base=$(basename "$r1_file")
    base=${base/_R1_*/}
    
    trimmed_r1_paired="$TRIMMED_DIR/${base}_R1_trimmed_paired.fastq.gz"
    trimmed_r1_unpaired="$TRIMMED_DIR/${base}_R1_trimmed_unpaired.fastq.gz"
    trimmed_r2_paired="$TRIMMED_DIR/${base}_R2_trimmed_paired.fastq.gz"
    trimmed_r2_unpaired="$TRIMMED_DIR/${base}_R2_trimmed_unpaired.fastq.gz"

    if [ ! -f "$trimmed_r1_paired" ]; then
        echo "Trimming $base..."
        if ! java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -threads "$THREADS" -phred33 \
            "$r1_file" "$r2_file" \
            "$trimmed_r1_paired" "$trimmed_r1_unpaired" \
            "$trimmed_r2_paired" "$trimmed_r2_unpaired" \
            TRAILING:10; then
            echo "ERROR: Trimmomatic failed for $base"
            rm -f "$trimmed_r1_paired" "$trimmed_r1_unpaired" \
                  "$trimmed_r2_paired" "$trimmed_r2_unpaired"
            exit 1
        fi

        # Validate output files
        for f in "$trimmed_r1_paired" "$trimmed_r2_paired"; do
            if ! gzip -t "$f"; then
                echo "ERROR: Trimmed file $f is corrupted"
                exit 1
            fi
        done
    fi
done

# STEP 2: Quality control
for file in "$TRIMMED_DIR"/*_trimmed_*.fastq.gz; do
    qc_report="$QC_DIR/$(basename "$file" .fastq.gz)_fastqc.html"
    if [ ! -f "$qc_report" ]; then
        echo "Running FastQC on $file..."
        fastqc -o "$QC_DIR" "$file"
    fi
done

# STEP 3: Paired-end alignment (corrected section)
index_files=("$HISAT2_INDEX".*.ht2)
if [ ! -f "${index_files[0]}" ]; then
    echo "ERROR: HISAT2 index not found at $HISAT2_INDEX"
    exit 1
fi

for r1_file in "$TRIMMED_DIR"/*_R1_trimmed_paired.fastq.gz; do
    base=$(basename "$r1_file" _R1_trimmed_paired.fastq.gz)
    r2_file="$TRIMMED_DIR/${base}_R2_trimmed_paired.fastq.gz"
    check_file_exists "$r2_file"
    
    bam="$ALIGN_DIR/${base}.bam"
    if [ ! -f "$bam" ]; then
        echo "Aligning $base..."
        strand=$(get_strandedness "$base" "HISAT2")
        hisat2_cmd="hisat2 -q -x $HISAT2_INDEX -1 $r1_file -2 $r2_file -p $THREADS"
        [ -n "$strand" ] && hisat2_cmd+=" --rna-strandness $strand"
        
        # Modified alignment pipeline with temp file and validation
        if ! $hisat2_cmd \
            | samtools view -@ $THREADS -bh \
            | samtools sort -@ $THREADS -o "$bam.tmp" - 
        then
            echo "ERROR: Alignment failed for $base" >&2
            rm -f "$bam.tmp"
            exit 1
        fi
        
        # Finalize BAM file
        mv "$bam.tmp" "$bam"
        
        # Validate BAM file
        samtools quickcheck "$bam" || { 
            echo "ERROR: Corrupted BAM $bam" >&2
            exit 1
        }
        samtools index -@ "$THREADS" "$bam"
    fi
done

# STEP 4: Quantification with featureCounts (corrected strandedness handling)
if [ ! -f "$QUANT_DIR/featurecounts_results.txt" ]; then
    echo "Running featureCounts..."
    samples=()
    strands=()

    while read -r sample strand; do
        bam_path="$ALIGN_DIR/${sample}.bam"
        check_file_exists "$bam_path"
        samples+=("$bam_path")
        strands+=("$(get_strandedness "$sample" "featureCounts")")
    done < "$STRANDEDNESS_FILE"

    # Get strandedness from first sample (assuming consistent across all)
    read -r first_sample first_strand <<< "$(head -n1 "$STRANDEDNESS_FILE")"
    fc_strand=$(get_strandedness "$first_sample" "featureCounts")

    featureCounts -T "$THREADS" -p -s "$fc_strand" -a "$ANNOTATION" \
        -o "$QUANT_DIR/featurecounts_results.txt" "${samples[@]}" || {
            echo "featureCounts failed"; exit 1;
        }
fi

# STEP 5: HTSeq counting
while read -r sample _; do
    bam="$ALIGN_DIR/${sample}.bam"
    output="$QUANT_DIR/${sample}_htseq.txt"
    
    check_file_exists "$bam"

    if [ ! -f "$output" ]; then
        echo "Running HTSeq on $sample..."
        strand=$(get_strandedness "$sample" "HTSeq")
        htseq-count -s "$strand" -f bam -r pos -t exon -i gene_id \
            "$bam" "$ANNOTATION" > "$output" || {
                echo "HTSeq failed for $sample"; exit 1;
            }
    fi
done < "$STRANDEDNESS_FILE"

# Final report
echo "
Pipeline complete!
QC reports:      $QC_DIR
Alignment files: $ALIGN_DIR
Quant results:   $QUANT_DIR
Total runtime:   $(($SECONDS / 60))m $(($SECONDS % 60))s
"