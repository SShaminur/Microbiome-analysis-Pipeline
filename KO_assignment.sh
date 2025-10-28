#!/bin/bash
# Unified KO Mapping and Processing Pipeline - Final Robust Version

# ======================
# CONFIGURATION
# ======================
set -euo pipefail  # Strict error handling

# Directory setup
BASE_DIR="/home/user/SR/SA-NGS/Chicken_data/"
INPUT_DIR="$BASE_DIR/RES/humann_results/raw_counts/"
OUTPUT_DIR="$BASE_DIR/F_KO"
TEMP_DIR="$BASE_DIR/temp_processing_$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$OUTPUT_DIR/processing_log_$(date +%Y%m%d_%H%M%S).txt"

# External resources
KO_MAPPING_FILE="/db/KO/map_ko_uniref90.txt"

# ======================
# INITIALIZATION
# ======================
initialize() {
    echo "Initializing KO processing pipeline..."
    
    # Check for bc command
    if ! command -v bc &> /dev/null; then
        echo "NOTE: 'bc' command not found, using awk for percentage calculations"
        USE_BC=false
    else
        USE_BC=true
    fi
    
    # Create directories
    mkdir -p "$OUTPUT_DIR" "$TEMP_DIR" || {
        echo "ERROR: Failed to create output directories" >&2
        exit 1
    }
    
    # Initialize log file
    exec > >(tee -a "$LOG_FILE") 2>&1
    echo "==== KO PROCESSING LOG - $(date) ===="
    echo "Input Directory: $INPUT_DIR"
    echo "Output Directory: $OUTPUT_DIR"
    echo "Temp Directory: $TEMP_DIR"
    echo "KO Mapping File: $KO_MAPPING_FILE"
    echo "===================================="
}

# ======================
# FUNCTION: CLEANUP
# ======================
cleanup() {
    local keep_temp=${1:-}
    
    if [[ -z "$keep_temp" && -d "$TEMP_DIR" ]]; then
        echo "Cleaning up temporary files..."
        rm -rf "$TEMP_DIR"
    else
        echo "Keeping temporary files in $TEMP_DIR"
    fi
}

# ======================
# FUNCTION: VALIDATE INPUTS
# ======================
validate_inputs() {
    echo "Validating input files..."
    
    # Check input directory exists
    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "ERROR: Input directory $INPUT_DIR not found" >&2
        exit 1
    fi
    
    # Check for at least one input file
    local input_files=("$INPUT_DIR"/*_raw_genefamilies.tsv)
    if [[ ${#input_files[@]} -eq 0 ]]; then
        echo "ERROR: No input files found in $INPUT_DIR" >&2
        exit 1
    fi
    
    # Check KO mapping file exists
    if [[ ! -f "$KO_MAPPING_FILE" ]]; then
        echo "ERROR: KO mapping file $KO_MAPPING_FILE not found" >&2
        exit 1
    fi
    
    echo "Found input files:"
    printf "  - %s\n" "${input_files[@]}"
    echo "Using KO mapping file: $KO_MAPPING_FILE"
}

# ======================
# FUNCTION: CALCULATE PERCENTAGE
# ======================
calculate_percentage() {
    local numerator=$1
    local denominator=$2
    
    if $USE_BC; then
        echo "scale=2; $numerator * 100 / $denominator" | bc
    else
        awk -v num="$numerator" -v den="$denominator" 'BEGIN {printf "%.2f", (num/den)*100}'
    fi
}

# ======================
# STEP 1: PROCESS GENE FAMILIES
# ======================
process_gene_families() {
    local input_file="$1"
    local sample_name=$(basename "$input_file" | cut -d'_' -f1)
    local output_file="$TEMP_DIR/${sample_name}_consolidated.tsv"
    
    echo "Processing $sample_name gene families..."
    
    # Consolidate gene families and remove taxonomy annotations
    awk -F'\t' '
    BEGIN {OFS="\t"}
    {
        # Process header
        if ($1 == "# Gene Family") {
            print $1, "'"${sample_name}"'_Abundance";
            next;
        }
        
        # Remove taxonomy annotations and sum abundances
        split($1, parts, "|");
        uniref_id = parts[1];
        abundance[uniref_id] += $2;
    }
    END {
        # Print consolidated results
        for (id in abundance) {
            # Format abundance to 5 decimal places
            printf "%s\t%.5f\n", id, abundance[id];
        }
    }' "$input_file" | sort > "$output_file" || {
        echo "ERROR: Failed to process gene families for $sample_name" >&2
        return 1
    }
    
    echo "Consolidated $(wc -l < "$output_file") unique gene families for $sample_name"
}

# ======================
# STEP 2: ASSIGN KO NUMBERS
# ======================
assign_ko_numbers() {
    local input_file="$1"
    local sample_name=$(basename "$input_file" | cut -d'_' -f1)
    local output_file="$TEMP_DIR/${sample_name}_with_KO.tsv"
    
    echo "Assigning KO numbers for $sample_name..."
    
    # First create optimized KO mapping file
    local ko_map_file="$TEMP_DIR/uniref_to_ko.tsv"
    awk '{
        ko = $1;
        for (i=2; i<=NF; i++) {
            print $i "\t" ko;
        }
    }' "$KO_MAPPING_FILE" > "$ko_map_file" || {
        echo "ERROR: Failed to process KO mapping file" >&2
        return 1
    }
    
    # Assign KO numbers to gene families
    awk -F'\t' '
    BEGIN {
        OFS="\t";
        # Load KO mappings
        while ((getline < "'"$ko_map_file"'") > 0) {
            if ($1 in ko_map) {
                ko_map[$1] = ko_map[$1] ";" $2;
            } else {
                ko_map[$1] = $2;
            }
        }
    }
    {
        # Process header
        if (NR == 1) {
            print $1, "KO_ID", $2;
            next;
        }
        
        # Print gene family with KO assignments
        print $1, ($1 in ko_map ? ko_map[$1] : ""), $2;
    }' "$input_file" > "$output_file" || {
        echo "ERROR: Failed to assign KO numbers" >&2
        return 1
    }
    
    # Calculate KO coverage statistics
    local total_genes=$(tail -n +2 "$output_file" | wc -l)
    local genes_with_ko=$(awk -F'\t' 'NR>1 && $2!=""' "$output_file" | wc -l)
    local ko_coverage=$(calculate_percentage "$genes_with_ko" "$total_genes")
    
    echo "KO assignment complete for $sample_name:"
    echo "- Total genes: $total_genes"
    echo "- Genes with KO: $genes_with_ko"
    echo "- KO coverage: $ko_coverage%"
}

# ======================
# STEP 3: PROCESS KO ABUNDANCES
# ======================
process_ko_abundances() {
    local input_file="$1"
    local sample_name=$(basename "$input_file" | cut -d'_' -f1)
    local output_file="$OUTPUT_DIR/${sample_name}_KO.tsv"
    
    echo "Processing KO abundances for $sample_name..."
    
    # Split multi-KO entries and consolidate abundances
    awk -F'\t' '
    BEGIN {OFS="\t"}
    NR == 1 {
        # Capture abundance column name
        abundance_col = $3;
        next;
    }
    $2 != "" {
        # Split multiple KO assignments
        split($2, kos, ";");
        for (i in kos) {
            ko = kos[i];
            sum[ko] += $3;
        }
    }
    END {
        print "KO_ID", "'"$sample_name"'";
        # Sort keys by KO ID (compatible with standard awk)
        count = 0;
        for (ko in sum) {
            keys[count++] = ko;
        }
        # Simple bubble sort for KO IDs
        for (i = 0; i < count; i++) {
            for (j = i+1; j < count; j++) {
                if (keys[i] > keys[j]) {
                    tmp = keys[i];
                    keys[i] = keys[j];
                    keys[j] = tmp;
                }
            }
        }
        # Output sorted results
        for (i = 0; i < count; i++) {
            printf "%s\t%.5f\n", keys[i], sum[keys[i]];
        }
    }' "$input_file" > "$output_file" || {
        echo "ERROR: Failed to process KO abundances" >&2
        return 1
    }
    
    echo "Created KO abundance file: $output_file"
}

# ======================
# STEP 4: MERGE SAMPLES
# ======================
merge_samples() {
    local samples=("$@")
    local output_file="$OUTPUT_DIR/merged_KO.tsv"
    
    echo "Merging ${#samples[@]} samples..."
    
    # Create a temporary directory for merge processing
    local merge_temp="$TEMP_DIR/merge_temp"
    mkdir -p "$merge_temp"
    
    # First, collect all unique KO IDs from all files
    awk -F'\t' '
    FNR == 1 {next}  # Skip headers
    {
        ko_ids[$1] = 1
    }
    END {
        for (ko in ko_ids) {
            print ko
        }
    }' "${samples[@]}" | sort -u > "$merge_temp/all_ko_ids.txt"
    
    # Initialize the output file with header
    printf "KO_ID" > "$output_file"
    for sample in "${samples[@]}"; do
        sample_name=$(basename "$sample" | cut -d'_' -f1)
        printf "\t%s" "$sample_name" >> "$output_file"
    done
    printf "\n" >> "$output_file"
    
    # Process each KO ID and merge abundances
    while read -r ko_id; do
        printf "%s" "$ko_id" >> "$output_file"
        for sample in "${samples[@]}"; do
            # Look up the abundance for this KO in each sample file
            abundance=$(awk -F'\t' -v ko="$ko_id" '$1 == ko {printf "%.5f", $2; exit}' "$sample" || echo "0")
            printf "\t%s" "$abundance" >> "$output_file"
        done
        printf "\n" >> "$output_file"
    done < "$merge_temp/all_ko_ids.txt"
    
    echo "Successfully created merged KO table: $output_file"
}

# ======================
# STEP 5: CREATE ZERO-FILLED TABLE
# ======================
create_zero_filled_table() {
    local input_file="$OUTPUT_DIR/merged_KO.tsv"
    local output_file="$OUTPUT_DIR/Zero_merged_KO.tsv"
    
    echo "Creating zero-filled version of merged KO table..."
    
    # Process the merged file to create zero-filled version
    awk -F'\t' '
    BEGIN {OFS="\t"}
    NR == 1 {
        # Print header without the KO_ID column
        for (i=2; i<=NF; i++) {
            printf "%s%s", $i, (i==NF ? "\n" : "\t");
        }
        next;
    }
    {
        # Process each data row
        for (i=2; i<=NF; i++) {
            # Replace empty values with 0, keep existing values
            printf "%s%s", ($i == "" || $i == "\t" ? "0" : $i), (i==NF ? "\n" : "\t");
        }
    }' "$input_file" > "$output_file" || {
        echo "ERROR: Failed to create zero-filled table" >&2
        return 1
    }
    
    echo "Created zero-filled KO table: $output_file"
}

# ======================
# MAIN EXECUTION
# ======================
main() {
    # Initialize environment
    initialize
    validate_inputs
    
    # Process each input file
    local processed_samples=()
    for input_file in "$INPUT_DIR"/*_raw_genefamilies.tsv; do
        local sample_name=$(basename "$input_file" | cut -d'_' -f1)
        
        echo "===================================="
        echo "Processing sample: $sample_name"
        echo "Input file: $input_file"
        
        # Step 1: Process gene families
        process_gene_families "$input_file" || {
            echo "Skipping $sample_name due to errors"
            continue
        }
        
        # Step 2: Assign KO numbers
        assign_ko_numbers "$TEMP_DIR/${sample_name}_consolidated.tsv" || {
            echo "Skipping $sample_name due to KO assignment errors"
            continue
        }
        
        # Step 3: Process KO abundances
        process_ko_abundances "$TEMP_DIR/${sample_name}_with_KO.tsv" && \
            processed_samples+=("$OUTPUT_DIR/${sample_name}_KO.tsv")
    done
    
    # Step 4: Merge samples if we have at least 2
    if [[ ${#processed_samples[@]} -ge 2 ]]; then
        merge_samples "${processed_samples[@]}" || {
            echo "Warning: Merge failed, keeping intermediate files"
            cleanup "keep"
            exit 1
        }
        
        # Step 5: Create zero-filled table
        create_zero_filled_table || {
            echo "Warning: Zero-filled table creation failed"
            cleanup "keep"
            exit 1
        }
    else
        echo "Note: Only ${#processed_samples[@]} samples processed - skipping merge and zero-filled table"
    fi
    
    # Final cleanup
    cleanup
    echo "===================================="
    echo "Processing complete! Results in $OUTPUT_DIR"
    echo "Log file: $LOG_FILE"
}

# Execute main function
main
