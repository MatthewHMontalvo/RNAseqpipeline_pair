#!/bin/bash

# Directory configuration
PROJECT_DIR="/Users/holdenmontalvo/RNAseq_pipeline"
RAW_DIR="$PROJECT_DIR/data/raw_fastq"
INDEX_DIR="$PROJECT_DIR/salmon_index"
OUTPUT_DIR="$PROJECT_DIR/data/salmon_output"
OUTPUT_FILE="$PROJECT_DIR/strandedness.txt"

# Create output directories
mkdir -p "$OUTPUT_DIR"

# Clear existing strandedness file
> "$OUTPUT_FILE"

# Process each R1 FASTQ pair
for R1_FILE in "$RAW_DIR"/*_R1_*.fastq.gz; do
    # Extract base sample name
    SAMPLE_BASE=$(basename "$R1_FILE" | sed 's/_R1_.*\.fastq\.gz//')
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    
    # Validate pair existence
    if [ ! -f "$R2_FILE" ]; then
        echo "ERROR: Missing pair for $R1_FILE"
        exit 1
    fi

    SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_BASE"
    mkdir -p "$SAMPLE_OUTPUT"
    
    echo "Processing $SAMPLE_BASE..."
    
    # Run Salmon with paired-end parameters
    salmon quant -i "$INDEX_DIR" \
        -l A \
        -1 "$R1_FILE" \
        -2 "$R2_FILE" \
        --validateMappings \
        --threads 10 \
        -o "$SAMPLE_OUTPUT"

    # Extract strandedness information
    STRANDEDNESS=$(grep '"expected_format"' "$SAMPLE_OUTPUT"/lib_format_counts.json | 
                   awk -F '"' '{print $4}' | 
                   sed -e 's/ISF/SF/' -e 's/ISR/SR/' -e 's/IU/U/')
    
    [ -z "$STRANDEDNESS" ] && STRANDEDNESS="U"

    # Save to strandedness file
    echo -e "$SAMPLE_BASE\t$STRANDEDNESS" >> "$OUTPUT_FILE"
    echo "Detected strandedness for $SAMPLE_BASE: $STRANDEDNESS"
done

echo "Strandedness detection complete. Results saved to $OUTPUT_FILE"