#!/bin/bash
# HUMAnN3 Pipeline for Multiple Samples - Parallel Processing with Cleanup
# Usage: ./humann_pipeline.sh <input_dir> <output_dir> <threads_per_sample> <max_parallel_samples>

set -euo pipefail
# Version: 1.6
# Author: Md Shaminur Rahman
# Date: Aug-2025

### PARAMETERS ###
INPUT_DIR="$1"
OUTPUT_DIR="$2"
THREADS_PER_SAMPLE="$3"
MAX_PARALLEL_SAMPLES="${4:-4}"  # Default to 4 parallel samples if not specified

### ALIGNMENT THRESHOLDS ###
IDENTITY_THRESHOLD="90.0"
QUERY_COVERAGE="90.0"
SUBJECT_COVERAGE="90.0"

### ADAPTERS ###
ADAPTERS="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"

### DATABASE PATHS ###
UNIREF90_DB="/db/uniref90/uniref"
CHOCOPHLAN_DB="/db/chocophlan/chocophlan/"
METAPHLAN_DB="/db/metaphlan4"
UTILITY_MAPPING="/db/uniref90/utility_mapping"

### PROGRESS TRACKING ###
PROGRESS_FILE="${OUTPUT_DIR}/progress.log"
touch "$PROGRESS_FILE"

### LOGGING SETUP ###
setup_logging() {
    LOG_DIR="${OUTPUT_DIR}/logs"
    mkdir -p "$LOG_DIR"
    exec > >(tee -a "${LOG_DIR}/pipeline.log") 2>&1
    echo "[$(date)] Logging initialized"
}

update_progress() {
    local SAMPLE=$1
    local STAGE=$2
    local STATUS=$3
    echo "$(date) | $SAMPLE | $STAGE | $STATUS" >> "$PROGRESS_FILE"
}

sample_log() {
    local SAMPLE=$1
    local MESSAGE="$2"
    echo "[$(date)] $MESSAGE" >> "${LOG_DIR}/${SAMPLE}.log"
}

### ERROR HANDLING ###
error_exit() {
    local LINE="$1"
    local MESSAGE="$2"
    echo "[$(date)] ERROR at line $LINE: $MESSAGE" >&2
    exit 1
}

trap 'error_exit $LINENO "Unexpected error occurred"' ERR

### TOOL VALIDATION ###
validate_tools() {
    echo "[$(date)] Validating required tools..."
    local REQUIRED_TOOLS=(
        trimmomatic
        repair.sh
        humann
        humann_join_tables
        bowtie2
        diamond
        parallel
    )
    
    local MISSING_TOOLS=()
    
    for tool in "${REQUIRED_TOOLS[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            MISSING_TOOLS+=("$tool")
        else
            echo "[$(date)] Found $tool: $(command -v "$tool")"
        fi
    done
    
    if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
        echo "[$(date)] ERROR: Missing required tools:" >&2
        printf "  - %s\n" "${MISSING_TOOLS[@]}" >&2
        echo "Install with:" >&2
        echo "  conda install -c biobakery humann metaphlan" >&2
        echo "  conda install -c bioconda bbmap bowtie2 diamond trimmomatic parallel" >&2
        exit 1
    fi
    
    # Verify adapter file exists
    if [ ! -f "$ADAPTERS" ]; then
        echo "[$(date)] ERROR: Adapter file not found: $ADAPTERS" >&2
        echo "Reinstall Trimmomatic with:" >&2
        echo "  conda install -c bioconda trimmomatic" >&2
        exit 1
    fi
    
    echo "[$(date)] All tools and adapters validated successfully"
}

### DATABASE VALIDATION ###
validate_databases() {
    local MISSING_DBS=()
    
    [ -d "$UNIREF90_DB" ] || MISSING_DBS+=("UniRef90")
    [ -d "$CHOCOPHLAN_DB" ] || MISSING_DBS+=("ChocoPhlAn")
    [ -d "$METAPHLAN_DB" ] || MISSING_DBS+=("MetaPhlAn")
    [ -d "$UTILITY_MAPPING" ] || MISSING_DBS+=("Utility Mapping")
    
    if [ ${#MISSING_DBS[@]} -gt 0 ]; then
        echo "[$(date)] ERROR: Missing required databases:" >&2
        printf "  - %s\n" "${MISSING_DBS[@]}" >&2
        echo "Download databases with:" >&2
        echo "  humann_databases --download uniref uniref90_diamond /db/uniref90" >&2
        echo "  humann_databases --download chocophlan full /db/chocophlan" >&2
        echo "  metaphlan --install --bowtie2db /db/metaphlan4 --index latest" >&2
        exit 1
    fi
    
    echo "[$(date)] All required databases are available"
}

### PARAMETER VALIDATION ###
validate_parameters() {
    [ $# -ge 3 ] || { 
        echo "Usage: $0 <input_dir> <output_dir> <threads_per_sample> [max_parallel_samples]" >&2
        echo "Note: Input directory should contain paired-end FASTQ files (R1 and R2)" >&2
        echo "      max_parallel_samples defaults to 4 if not specified" >&2
        exit 1
    }
    
    [ -d "$INPUT_DIR" ] || { echo "[$(date)] ERROR: Input directory not found: $INPUT_DIR" >&2; exit 1; }
    
    [[ "$THREADS_PER_SAMPLE" =~ ^[0-9]+$ ]] && [ "$THREADS_PER_SAMPLE" -gt 0 ] || { 
        echo "[$(date)] ERROR: Invalid thread count per sample: $THREADS_PER_SAMPLE" >&2; exit 1
    }
    
    [[ "$MAX_PARALLEL_SAMPLES" =~ ^[0-9]+$ ]] && [ "$MAX_PARALLEL_SAMPLES" -gt 0 ] || { 
        echo "[$(date)] ERROR: Invalid max parallel samples: $MAX_PARALLEL_SAMPLES" >&2; exit 1
    }
    
    # Calculate total threads needed
    local TOTAL_THREADS=$((THREADS_PER_SAMPLE * MAX_PARALLEL_SAMPLES))
    echo "[$(date)] Configuration:"
    echo "[$(date)]   Threads per sample: $THREADS_PER_SAMPLE"
    echo "[$(date)]   Max parallel samples: $MAX_PARALLEL_SAMPLES"
    echo "[$(date)]   Total threads required: $TOTAL_THREADS"
}

### INITIALIZATION ###
initialize() {
    mkdir -p "$OUTPUT_DIR" || { echo "[$(date)] ERROR: Cannot create output directory" >&2; exit 1; }
    mkdir -p "${OUTPUT_DIR}/"{trimmed,interleaved,humann_results/raw_counts,logs,tmp} || { 
        echo "[$(date)] ERROR: Cannot create subdirectories" >&2; exit 1 
    }
    
    echo "[$(date)] ===== STARTING HUMAnN3 PIPELINE ====="
    echo "[$(date)] Parameters:"
    echo "[$(date)]   Input Directory: $INPUT_DIR"
    echo "[$(date)]   Output Directory: $OUTPUT_DIR"
    echo "[$(date)]   Threads per Sample: $THREADS_PER_SAMPLE"
    echo "[$(date)]   Max Parallel Samples: $MAX_PARALLEL_SAMPLES"
    echo "[$(date)] Alignment Thresholds:"
    echo "[$(date)]   Identity: $IDENTITY_THRESHOLD%"
    echo "[$(date)]   Query Coverage: $QUERY_COVERAGE%"
    echo "[$(date)]   Subject Coverage: $SUBJECT_COVERAGE%"
}

### GET SAMPLE PAIRS ###
get_sample_pairs() {
    # Find all R1 files in the input directory
    for r1 in "${INPUT_DIR}"/*_R1*.f*q*; do
        # Skip if not a file (handles case where no files match)
        [ -f "$r1" ] || continue
        
        # Extract base sample name (remove _R1 and everything after)
        sample_name=$(basename "$r1" | sed -E 's/_R1.*//')
        
        # Find matching R2 file - try multiple patterns
        r2=""
        # Try exact _R1 to _R2 replacement first
        possible_r2="${r1/_R1/_R2}"
        if [ -f "$possible_r2" ]; then
            r2="$possible_r2"
        else
            # Try alternative patterns if exact replacement doesn't work
            possible_r2=$(echo "$r1" | sed -E 's/_R1([._])/_R2\1/')
            if [ -f "$possible_r2" ]; then
                r2="$possible_r2"
            else
                echo "[$(date)] ERROR: Could not find R2 file for $r1" >&2
                continue
            fi
        fi
        
        echo "$sample_name:$r1:$r2"
    done | sort -u
}

### RUN HUMANN ###
run_humann() {
    local INPUT="$1"
    local SAMPLE="$2"
    
    echo "=== Currently processing sample: $SAMPLE ==="
    update_progress "$SAMPLE" "HUMAnN" "STARTED"
    
    humann \
        --input "$INPUT" \
        --output "${OUTPUT_DIR}/humann_results" \
        --threads "$THREADS_PER_SAMPLE" \
        --protein-database "$UNIREF90_DB" \
        --nucleotide-database "$CHOCOPHLAN_DB" \
        --metaphlan-options "--bowtie2db $METAPHLAN_DB" \
        --translated-alignment diamond \
        --translated-identity-threshold "$IDENTITY_THRESHOLD" \
        --translated-query-coverage-threshold "$QUERY_COVERAGE" \
        --translated-subject-coverage-threshold "$SUBJECT_COVERAGE" \
        --search-mode uniref90 \
        --pathways metacyc \
        --gap-fill on \
        --minpath on \
        --output-basename "$SAMPLE" \
        --output-format tsv \
        --memory-use maximum \
        --verbose
    
    if [ $? -eq 0 ]; then
        update_progress "$SAMPLE" "HUMAnN" "COMPLETED"
    else
        update_progress "$SAMPLE" "HUMAnN" "FAILED"
        return 1
    fi
}

### PROCESS SAMPLE ###
process_sample() {
    local R1="$1"
    local R2="$2"
    local SAMPLE="$3"
    local TMP_DIR="${OUTPUT_DIR}/tmp/${SAMPLE}"
    
    # Create sample-specific directories and log file immediately
    mkdir -p "$TMP_DIR"
    exec > >(tee -a "${LOG_DIR}/${SAMPLE}.log") 2>&1
    
    echo "[$(date)] Starting processing for sample: $SAMPLE"
    update_progress "$SAMPLE" "INIT" "STARTED"
    echo "[$(date)] Creating temporary directory: $TMP_DIR"

    ### 1. QUALITY TRIMMING ###
    echo "[$(date)] Running Trimmomatic for $SAMPLE"
    update_progress "$SAMPLE" "TRIMMING" "STARTED"
    trimmomatic PE -threads "$THREADS_PER_SAMPLE" \
        "$R1" "$R2" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_paired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_unpaired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_paired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_unpaired.fq.gz" \
        "ILLUMINACLIP:${ADAPTERS}:2:30:10:3:TRUE" \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36 || {
        echo "[$(date)] ERROR: Trimmomatic failed for $SAMPLE" >&2
        update_progress "$SAMPLE" "TRIMMING" "FAILED"
        return 1
    }
    update_progress "$SAMPLE" "TRIMMING" "COMPLETED"
    
    ### 2. CREATE INTERLEAVED FASTQ ###
    echo "[$(date)] Running repair.sh (BBTools) for $SAMPLE"
    update_progress "$SAMPLE" "INTERLEAVING" "STARTED"
    repair.sh \
        in1="${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_paired.fq.gz" \
        in2="${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_paired.fq.gz" \
        out="${OUTPUT_DIR}/interleaved/${SAMPLE}_interleaved.fq.gz" || {
        echo "[$(date)] ERROR: Interleaving failed for $SAMPLE" >&2
        update_progress "$SAMPLE" "INTERLEAVING" "FAILED"
        return 1
    }
    update_progress "$SAMPLE" "INTERLEAVING" "COMPLETED"
    
    ### 3. RUN HUMAnN3 ###
    echo "[$(date)] Running HUMAnN3 for $SAMPLE"
    run_humann "${OUTPUT_DIR}/interleaved/${SAMPLE}_interleaved.fq.gz" "$SAMPLE" || {
        echo "[$(date)] ERROR: HUMAnN3 failed for $SAMPLE" >&2
        return 1
    }
    
    ### 4. PRESERVE RAW COUNTS ###
    echo "[$(date)] Saving raw counts for $SAMPLE"
    update_progress "$SAMPLE" "RAW_COUNTS" "STARTED"
    for file in "${OUTPUT_DIR}/humann_results/${SAMPLE}"*genefamilies.tsv; do
        if [ -f "$file" ]; then
            cp "$file" "${OUTPUT_DIR}/humann_results/raw_counts/${SAMPLE}_raw_genefamilies.tsv"
        else
            echo "[$(date)] WARNING: No genefamilies file found"
            update_progress "$SAMPLE" "RAW_COUNTS" "WARNING: No genefamilies"
        fi
    done
    update_progress "$SAMPLE" "RAW_COUNTS" "COMPLETED"
    
    ### 5. CLEAN UP INTERMEDIATE FILES ###
    echo "[$(date)] Cleaning up intermediate files for $SAMPLE"
    update_progress "$SAMPLE" "CLEANUP" "STARTED"
    rm -f \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_paired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R1_unpaired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_paired.fq.gz" \
        "${OUTPUT_DIR}/trimmed/${SAMPLE}_R2_unpaired.fq.gz" \
        "${OUTPUT_DIR}/interleaved/${SAMPLE}_interleaved.fq.gz"
    update_progress "$SAMPLE" "CLEANUP" "COMPLETED"
    
    echo "[$(date)] Completed processing for sample: $SAMPLE"
    update_progress "$SAMPLE" "FINAL" "COMPLETED"
}

### CLEANUP HUMANN RESULTS ###
cleanup_humann_results() {
    echo "[$(date)] Cleaning up HUMAnN results directory..."
    update_progress "GLOBAL" "CLEANUP" "STARTED"
    
    # Keep only essential files
    find "${OUTPUT_DIR}/humann_results" -type f \( -name "*.tsv" -o -name "*.log" \) -print0 | while IFS= read -r -d '' file; do
        # Skip files in raw_counts directory
        if [[ "$file" != *"/raw_counts/"* ]]; then
            rm -f "$file"
        fi
    done
    
    # Remove empty directories except raw_counts and logs
    find "${OUTPUT_DIR}/humann_results" -mindepth 1 -type d ! -name "raw_counts" ! -name "logs" -exec rmdir --ignore-fail-on-non-empty {} \;
    
    echo "[$(date)] HUMAnN results cleanup completed"
    update_progress "GLOBAL" "CLEANUP" "COMPLETED"
}

### EXPORT FUNCTIONS FOR PARALLEL ###
export -f process_sample run_humann sample_log update_progress
export INPUT_DIR OUTPUT_DIR THREADS_PER_SAMPLE IDENTITY_THRESHOLD QUERY_COVERAGE SUBJECT_COVERAGE
export ADAPTERS UNIREF90_DB CHOCOPHLAN_DB METAPHLAN_DB UTILITY_MAPPING LOG_DIR PROGRESS_FILE

### MAIN PROCESSING ###
main() {
    setup_logging
    validate_parameters "$@"
    validate_tools
    validate_databases
    initialize
    
    # Process all samples
    local PROCESSED=0
    local FAILED=0
    
    # Get sample pairs
    echo "[$(date)] Searching for sample pairs in $INPUT_DIR"
    SAMPLE_PAIRS=()
    while IFS= read -r line; do
        [ -n "$line" ] && SAMPLE_PAIRS+=("$line")
    done < <(get_sample_pairs)
    
    # Check if any samples were found
    if [ ${#SAMPLE_PAIRS[@]} -eq 0 ]; then
        echo "[$(date)] ERROR: No valid sample pairs found in $INPUT_DIR" >&2
        echo "[$(date)] Expected files with patterns like:" >&2
        echo "[$(date)]   sample_R1.fastq.gz and sample_R2.fastq.gz" >&2
        echo "[$(date)]   sample_R1_001.fq and sample_R2_001.fq" >&2
        echo "[$(date)]   sample.R1.fq and sample.R2.fq" >&2
        exit 1
    fi
    
    echo "[$(date)] Found ${#SAMPLE_PAIRS[@]} sample pairs to process"
    
    # Process samples in parallel
    printf "%s\n" "${SAMPLE_PAIRS[@]}" | \
        parallel -j "$MAX_PARALLEL_SAMPLES" --colsep ':' --halt soon,fail=1 --joblog "${LOG_DIR}/parallel.log" \
        "process_sample {2} {3} {1}" || {
        echo "[$(date)] ERROR: Parallel processing failed for one or more samples" >&2
        FAILED=1
    }
    
    # Count processed samples
    PROCESSED=$(find "${OUTPUT_DIR}/humann_results" -maxdepth 1 -name "*genefamilies.tsv" | wc -l)
    
    ### CLEANUP HUMANN RESULTS ###
    cleanup_humann_results
    
    # Final status report
    echo "[$(date)] ===== PIPELINE COMPLETED ====="
    echo "[$(date)] Samples processed successfully: $PROCESSED"
    if [ $FAILED -ne 0 ]; then
        echo "[$(date)] ERROR: Some samples failed - check individual log files" >&2
    fi
    
    # Cleanup temporary files
    rm -rf "${OUTPUT_DIR}/tmp" 2>/dev/null
    
    if [ $FAILED -ne 0 ]; then
        exit 1
    fi
    
    # Launch KO processing pipeline if successful
    echo "[$(date)] Launching KO processing pipeline..."
    /bin/bash ./KO_assignment.sh "${OUTPUT_DIR}/humann_results/raw_counts/"
}

### EXECUTE ###
main "$@"
exit 0
