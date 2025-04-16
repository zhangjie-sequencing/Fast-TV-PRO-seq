#!/bin/bash

# Check required tools
command -v cutadapt >/dev/null 2>&1 || { echo >&2 "Error: cutadapt not installed. Aborting."; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "Error: bowtie2 not installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "Error: samtools not installed. Aborting."; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo >&2 "Error: bedtools not installed. Aborting."; exit 1; }

# ========== CONFIGURATION ==========
BOWTIE2_INDEX="/mnt/c/data/Bowtie2_WBcel235/WBcel235"  # ← MUST UPDATE THIS PATH!
# ===================================

# Processing function
process_fastq() {
    input_file="$1"
    # Fix filename handling (preserve original base name)
    base_name=$(basename "$input_file" | sed 's/\.f\(ast\)\?q\.\?gz$//')
    
    # Step 1: Adapter trimming
    echo "Trimming: $input_file"
    cutadapt -m 20 -e 0.05 --cores 20 -a TGGAATTCTCGGGTGCCAAGG \
        -o "${base_name}_trimmed.fq.gz" "$input_file" || exit 1
    
    # Step 2: Alignment
    echo "Aligning: $base_name"
    bowtie2 -p 20 -k 1 --very-sensitive --no-unal -x "$BOWTIE2_INDEX" -U "${base_name}_trimmed.fq.gz" -S "${base_name}.sam" || exit 1
    
    # Step 3: SAM → sorted BAM
    samtools view -@ 20 -bS "${base_name}.sam" > "${base_name}.bam" || exit 1
    samtools sort -@ 20 "${base_name}.bam" -o "${base_name}.sorted.bam" || exit 1
    
    # Step 4: Generate BedGraph
    bedtools genomecov -strand - -5 -bga -ibam "${base_name}.sorted.bam" > "${base_name}_p.bedgraph"
    bedtools genomecov -strand + -5 -bga -ibam "${base_name}.sorted.bam" > "${base_name}_m.bedgraph"
    
    # Cleanup with error suppression
    rm -f "${base_name}.sam" "${base_name}.bam" "${base_name}_trimmed.fq.gz"
    echo "Successfully processed: ${base_name}_[pm].bedgraph"
}

export -f process_fastq

# Handle interrupt signal
trap 'echo -e "\nInterrupted! Partial files may remain."; exit 130' SIGINT

# Find and process files
find . -type f \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) -print0 \
| parallel -0 -j 40 process_fastq

echo "All files processed successfully."