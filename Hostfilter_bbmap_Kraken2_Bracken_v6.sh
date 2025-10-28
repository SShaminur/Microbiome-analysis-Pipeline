#!/bin/bash

###############################################################################
# Kraken2/Bracken Analysis Pipeline with Trimmomatic Preprocessing and Host Filtering
# 
# Description:
#   This script performs metagenomic analysis using:
#   1. Trimmomatic for quality control
#   2. BBduk for host read filtering (optional)
#   3. Kraken2 for taxonomic classification
#   4. Bracken for abundance estimation
#   It also generates summary tables of the results.
#
# Usage:
#   ./Hostfilter_bbmap_Kraken2_Bracken_v6.sh [options] <input_files>
#
# Options:
#   -h, --help      Show this help message and exit
#   --skip-trim     Skip the trimming step (use existing files in paired/)
#   *.fastq.gz      Process files with the given pattern (must be paired-end)
#
# Input Requirements:
#   - Paired-end FASTQ files named with _R1.fastq.gz and _R2.fastq.gz suffixes
#
# Output:
#   - trimmed_results/: Trimmomatic logs
#   - paired/: Trimmed paired-end files
#   - host_filtered/: Host-filtered files and reports
#   - kraken2_results/raw/: Kraken2 classification results
#   - kraken2_results/bracken/: Bracken abundance estimates
#   - kraken2_results/Final_Species.tsv: Species-level summary
#   - kraken2_results/Final_Genus.tsv: Genus-level summary
#
# Dependencies:
#   - Trimmomatic, BBtools (bbduk), Kraken2, Bracken
#   - Properly configured Kraken2 database
#   - Host reference genome (if host filtering enabled)
#
# Version: 1.6
# Author: Md Shaminur Rahman
# Date: Aug-2025
###############################################################################

##################################################
### CONFIGURATION - EDIT THESE PARAMETERS FIRST ##
##################################################

# System resources
TOTAL_THREADS=80
THREADS_PER_TRIMMING_JOB=20
MEMORY_PER_THREAD="10G"

# Trimmomatic parameters
ADAPTERS_FILE="/home/user/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa"
HEADCROP=7
LEADING=20
TRAILING=20
SLIDINGWINDOW="20:20"
MINLEN=40

# Host filtering parameters
HOST_FILTER=true  # Set to false to disable host filtering
HOST_REFERENCE="/home/user/SR/db/chicken/ncbi_dataset/ncbi_dataset/data/GCA_000002315.5/GCA_000002315.5_GRCg6a_genomic.fna"

# Kraken2 parameters
KRAKEN_DB="/home/user/SR/db/kraken2/kraken_db/"
KRAKEN_CONFIDENCE=0.1
KRAKEN_MIN_HIT_GROUPS=2
KRAKEN_MIN_QUALITY=20

# Bracken parameters
BRACKEN_READ_LENGTH=150
BRACKEN_LEVEL="S"  # S for species, G for genus

##################################################
### MAIN SCRIPT - EDIT BELOW THIS LINE CAREFULLY #
##################################################

# Initialize color variables
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
reset=$(tput sgr0)

# Calculate maximum parallel jobs
MAX_PARALLEL_JOBS=$((TOTAL_THREADS / THREADS_PER_TRIMMING_JOB))
[ $MAX_PARALLEL_JOBS -lt 1 ] && MAX_PARALLEL_JOBS=1

usage() {
  echo -e "${green}
$(basename "$0") - Kraken2/Bracken Analysis Pipeline with Host Filtering

${yellow}Usage:${reset}
  $(basename "$0") [options] <input_files>

${yellow}Options:${reset}
  -h, --help      Show this help message and exit
  --skip-trim     Skip the trimming step (use existing files in paired/)
  *.fastq.gz      Process files with the given pattern (must are paired-end)

${yellow}Configuration:${reset}
  Threads: $TOTAL_THREADS total, $THREADS_PER_TRIMMING_JOB per trimming job
  Host filtering: $HOST_FILTER ${HOST_FILTER:+($HOST_REFERENCE)}
  Kraken DB: $KRAKEN_DB

${yellow}Output Folders:${reset}
  trimmed_results/ - Trimmomatic logs
  paired/ - Trimmed paired-end files
  host_filtered/ - Host-filtered files and reports
  kraken2_results/ - Kraken2 and Bracken results
  kraken2_results/Final_Species.tsv - Species-level summary
  kraken2_results/Final_Genus.tsv - Genus-level summary
${reset}"
  exit 0
}

create_directories() {
  echo -e "${yellow}Creating required directories...${reset}"
  mkdir -p trimmed_results paired host_filtered/status kraken2_results/{raw,bracken/{species,genus}}
}

count_reads() {
  local file=$1
  if [[ $file == *.gz ]]; then
    zcat "$file" | awk 'END {print NR/4}'
  else
    awk 'END {print NR/4}' "$file"
  fi
}

generate_summary_tables() {
  echo -e "${yellow}Generating summary tables...${reset}"
  
  # Process species files
  local species_files=(kraken2_results/bracken/species/*.tsv)
  if [ ${#species_files[@]} -gt 0 ]; then
    echo -e "${yellow}Processing species files...${reset}"
    
    # First pass to collect all taxonomy IDs and sample names
    declare -A tax_ids
    declare -A samples
    for file in "${species_files[@]}"; do
      sample=$(basename "$file" .tsv)
      samples["$sample"]=1
      
      while IFS=$'\t' read -r name tax_id level kraken_reads added_reads new_est_reads fraction; do
        if [[ "$tax_id" =~ ^[0-9]+$ ]]; then
          tax_ids["$tax_id"]=1
        fi
      done < <(tail -n +2 "$file")
    done
    
    # Generate header
    echo -n "taxonomy_id" > kraken2_results/Final_Species.tsv
    for sample in "${!samples[@]}"; do
      echo -ne "\t${sample}_kraken_assigned_reads" >> kraken2_results/Final_Species.tsv
    done
    echo >> kraken2_results/Final_Species.tsv
    
    # Generate data rows
    for tax_id in "${!tax_ids[@]}"; do
      echo -n "$tax_id" >> kraken2_results/Final_Species.tsv
      for sample in "${!samples[@]}"; do
        file="kraken2_results/bracken/species/${sample}.tsv"
        count=$(awk -v id="$tax_id" -F'\t' '$2 == id {print $4}' "$file" || echo "0")
        echo -ne "\t${count:-0}" >> kraken2_results/Final_Species.tsv
      done
      echo >> kraken2_results/Final_Species.tsv
    done
    
    echo -e "${green}Created species summary: kraken2_results/Final_Species.tsv${reset}"
  else
    echo -e "${red}No species files found to generate summary${reset}"
  fi
  
  # Process genus files (same approach as species)
  local genus_files=(kraken2_results/bracken/genus/*.tsv)
  if [ ${#genus_files[@]} -gt 0 ]; then
    echo -e "${yellow}Processing genus files...${reset}"
    
    # First pass to collect all taxonomy IDs and sample names
    declare -A tax_ids
    declare -A samples
    for file in "${genus_files[@]}"; do
      sample=$(basename "$file" .tsv)
      samples["$sample"]=1
      
      while IFS=$'\t' read -r name tax_id level kraken_reads added_reads new_est_reads fraction; do
        if [[ "$tax_id" =~ ^[0-9]+$ ]]; then
          tax_ids["$tax_id"]=1
        fi
      done < <(tail -n +2 "$file")
    done
    
    # Generate header
    echo -n "taxonomy_id" > kraken2_results/Final_Genus.tsv
    for sample in "${!samples[@]}"; do
      echo -ne "\t${sample}_kraken_assigned_reads" >> kraken2_results/Final_Genus.tsv
    done
    echo >> kraken2_results/Final_Genus.tsv
    
    # Generate data rows
    for tax_id in "${!tax_ids[@]}"; do
      echo -n "$tax_id" >> kraken2_results/Final_Genus.tsv
      for sample in "${!samples[@]}"; do
        file="kraken2_results/bracken/genus/${sample}.tsv"
        count=$(awk -v id="$tax_id" -F'\t' '$2 == id {print $4}' "$file" || echo "0")
        echo -ne "\t${count:-0}" >> kraken2_results/Final_Genus.tsv
      done
      echo >> kraken2_results/Final_Genus.tsv
    done
    
    echo -e "${green}Created genus summary: kraken2_results/Final_Genus.tsv${reset}"
  else
    echo -e "${red}No genus files found to generate summary${reset}"
  fi
}

run_trimmomatic() {
  local R1=$1
  local R2=$2
  local sample_name=$3
  
  R1_pair="paired/${sample_name}_R1_paired.fastq.gz"
  R1_unpair="paired/${sample_name}_R1_unpaired.fastq.gz"
  R2_pair="paired/${sample_name}_R2_paired.fastq.gz"
  R2_unpair="paired/${sample_name}_R2_unpaired.fastq.gz"
  
  echo -e "${yellow}[Trimming] ${sample_name}${reset}"
  
  trimmomatic PE \
    -threads $THREADS_PER_TRIMMING_JOB \
    -phred33 "$R1" "$R2" "$R1_pair" "$R1_unpair" "$R2_pair" "$R2_unpair" \
    HEADCROP:$HEADCROP \
    ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10:3:TRUE \
    LEADING:$LEADING \
    TRAILING:$TRAILING \
    SLIDINGWINDOW:$SLIDINGWINDOW \
    MINLEN:$MINLEN \
    &>> "trimmed_results/${sample_name}.log"
  
  # Remove unpaired files
  rm -f "$R1_unpair" "$R2_unpair"
  
  if [ $? -eq 0 ]; then
    echo -e "${green}Trimming completed for ${sample_name}${reset}"
  else
    echo -e "${red}Error trimming ${sample_name}${reset}"
    return 1
  fi
}

run_host_filtering() {
  local sample_name=$1
  local skip_trim=$2
  
  echo -e "${yellow}[Host Filtering] ${sample_name}${reset}"
  
  # Determine input files based on whether trimming was skipped
  if [ "$skip_trim" = true ]; then
    # Use original input files
    R1="${sample_name}_R1.fastq.gz"
    R2="${sample_name}_R2.fastq.gz"
  else
    # Use trimmed files
    R1="paired/${sample_name}_R1_paired.fastq.gz"
    R2="paired/${sample_name}_R2_paired.fastq.gz"
  fi
  
  # Check if input files exist
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo -e "${red}Error: Input files not found for ${sample_name}${reset}"
    echo -e "${red}Looking for: $R1 and $R2${reset}"
    return 1
  fi
  
  # Count input reads
  input_reads_R1=$(count_reads "$R1")
  input_reads_R2=$(count_reads "$R2")
  total_input_reads=$((input_reads_R1 + input_reads_R2))
  
  # Run BBduk for host filtering with optimized memory settings
  echo -e "${yellow}Running BBduk for host filtering...${reset}"
  
  # Calculate memory per BBduk job (80GB per job for 1TB total)
  BBDUK_MEMORY="80g"
  BBDUK_THREADS=$((TOTAL_THREADS / MAX_PARALLEL_JOBS))
  
  bbduk.sh -Xmx${BBDUK_MEMORY} \
    in1="$R1" \
    in2="$R2" \
    out1="host_filtered/${sample_name}_R1_filtered.fastq.gz" \
    out2="host_filtered/${sample_name}_R2_filtered.fastq.gz" \
    ref="$HOST_REFERENCE" \
    k=31 \
    threads=${BBDUK_THREADS} \
    stats="host_filtered/status/${sample_name}_bbduk_stats.txt" \
    prealloc \
    -da \
    &>> "host_filtered/status/${sample_name}_bbduk.log"
  
  if [ $? -ne 0 ]; then
    echo -e "${red}BBduk host filtering failed for ${sample_name}${reset}"
    echo -e "${yellow}Check log file: host_filtered/status/${sample_name}_bbduk.log${reset}"
    return 1
  fi
  
  # Count output reads
  output_reads_R1=$(count_reads "host_filtered/${sample_name}_R1_filtered.fastq.gz")
  output_reads_R2=$(count_reads "host_filtered/${sample_name}_R2_filtered.fastq.gz")
  total_output_reads=$((output_reads_R1 + output_reads_R2))
  
  # Calculate percentages
  filtered_reads=$((total_input_reads - total_output_reads))
  filtered_percent=$(awk -v f="$filtered_reads" -v t="$total_input_reads" 'BEGIN {printf "%.2f", (f/t)*100}')
  remaining_percent=$(awk -v o="$total_output_reads" -v t="$total_input_reads" 'BEGIN {printf "%.2f", (o/t)*100}')
  
  # Write status report
  {
    echo "Host Filtering Status for ${sample_name}"
    echo "====================================="
    echo "Input Reads: $total_input_reads"
    echo "Host Reads Filtered: $filtered_reads ($filtered_percent%)"
    echo "Unmapped Reads: $total_output_reads ($remaining_percent%)"
    echo "Filtered Files:"
    echo "- R1: host_filtered/${sample_name}_R1_filtered.fastq.gz"
    echo "- R2: host_filtered/${sample_name}_R2_filtered.fastq.gz"
    echo "====================================="
    echo "Generated on: $(date)"
  } > "host_filtered/status/${sample_name}_filtering_report.txt"
  
  echo -e "${green}Host filtering completed for ${sample_name}${reset}"
  echo -e "${yellow}Status report: host_filtered/status/${sample_name}_filtering_report.txt${reset}"
}

run_kraken2() {
  local sample_name=$1
  local skip_trim=$2
  
  echo -e "${yellow}[Kraken2] ${sample_name}${reset}"
  
  # Create sample-specific output directory
  mkdir -p "kraken2_results/raw/${sample_name}"
  
  # Determine input files based on whether host filtering was done
  if [ "$HOST_FILTER" = true ]; then
    R1="host_filtered/${sample_name}_R1_filtered.fastq.gz"
    R2="host_filtered/${sample_name}_R2_filtered.fastq.gz"
  elif [ "$skip_trim" = true ]; then
    # Use original input files
    R1="${sample_name}_R1.fastq.gz"
    R2="${sample_name}_R2.fastq.gz"
  else
    # Use trimmed files
    R1="paired/${sample_name}_R1_paired.fastq.gz"
    R2="paired/${sample_name}_R2_paired.fastq.gz"
  fi
  
  # Check if input files exist
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo -e "${red}Error: Input files not found for ${sample_name}${reset}"
    echo -e "${red}Looking for: $R1 and $R2${reset}"
    return 1
  fi
  
  local kraken_cmd=(
    kraken2
    --db "$KRAKEN_DB"
    --threads $((TOTAL_THREADS / MAX_PARALLEL_JOBS))
    --paired
    --confidence $KRAKEN_CONFIDENCE
    --minimum-hit-groups $KRAKEN_MIN_HIT_GROUPS
    --minimum-base-quality $KRAKEN_MIN_QUALITY
    --report "kraken2_results/raw/${sample_name}/report.tsv"
    --output "kraken2_results/raw/${sample_name}/classifications.txt"
    --use-names
    --gzip-compressed
    "$R1"
    "$R2"
  )
  
  echo -e "${yellow}Running: ${kraken_cmd[*]}${reset}"
  if ! "${kraken_cmd[@]}" &>> "kraken2_results/raw/${sample_name}/kraken.log"; then
    echo -e "${red}Error running Kraken2 on ${sample_name}${reset}"
    echo -e "${yellow}Check log file: kraken2_results/raw/${sample_name}/kraken.log${reset}"
    return 1
  fi
  
  echo -e "${green}Kraken2 completed for ${sample_name}${reset}"
}

run_bracken() {
  local sample_name=$1
  
  echo -e "${yellow}[Bracken] ${sample_name}${reset}"
  
  # Check if Kraken2 report exists
  if [ ! -f "kraken2_results/raw/${sample_name}/report.tsv" ]; then
    echo -e "${red}Error: Kraken2 report not found for ${sample_name}${reset}"
    return 1
  fi
  
  # Species level
  bracken \
    -d "$KRAKEN_DB" \
    -i "kraken2_results/raw/${sample_name}/report.tsv" \
    -o "kraken2_results/bracken/species/${sample_name}.tsv" \
    -l S \
    -r $BRACKEN_READ_LENGTH \
    -t $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) \
    &>> "kraken2_results/raw/${sample_name}/bracken.log"
  
  # Genus level
  bracken \
    -d "$KRAKEN_DB" \
    -i "kraken2_results/raw/${sample_name}/report.tsv" \
    -o "kraken2_results/bracken/genus/${sample_name}.tsv" \
    -l G \
    -r $BRACKEN_READ_LENGTH \
    -t $((TOTAL_THREADS / MAX_PARALLEL_JOBS)) \
    &>> "kraken2_results/raw/${sample_name}/bracken.log"
  
  if [ $? -eq 0 ]; then
    echo -e "${green}Bracken completed for ${sample_name}${reset}"
  else
    echo -e "${red}Error running Bracken on ${sample_name}${reset}"
    echo -e "${yellow}Check log file: kraken2_results/raw/${sample_name}/bracken.log${reset}"
    return 1
  fi
}

process_sample() {
  local sample_name=$1
  local skip_trim=$2
  
  # Step 1: Trimming (unless skipped)
  if [ "$skip_trim" = false ]; then
    run_trimmomatic "${sample_name}_R1.fastq.gz" "${sample_name}_R2.fastq.gz" "$sample_name" || return 1
  fi
  
  # Step 2: Host filtering (if enabled)
  if [ "$HOST_FILTER" = true ]; then
    run_host_filtering "$sample_name" "$skip_trim" || return 1
  fi
  
  # Step 3: Kraken2
  run_kraken2 "$sample_name" "$skip_trim" || return 1
  
  # Step 4: Bracken
  run_bracken "$sample_name" || return 1
  
  echo -e "${green}Completed analysis for ${sample_name}${reset}"
}

main() {
  local skip_trim=false
  
  # Parse arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help) usage ;;
      --skip-trim) 
        skip_trim=true
        shift ;;
      *)
        break ;;
    esac
  done

  if [ $# -eq 0 ]; then
    echo -e "${red}Error: No input files specified${reset}"
    usage
    exit 1
  fi

  create_directories
  
  # Process samples
  for r1_file in "$@"; do
    if [[ "$r1_file" == *_R1.fastq.gz ]]; then
      sample_name=$(basename "$r1_file" | sed 's/_R1\.fastq\.gz//')
      r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
      
      if [ ! -f "$r2_file" ]; then
        echo -e "${red}Error: Missing R2 file for sample ${sample_name}${reset}"
        continue
      fi
      
      echo -e "\n${green}Processing sample: ${sample_name}${reset}"
      process_sample "$sample_name" "$skip_trim" &
      
      # Limit number of parallel jobs
      if [[ $(jobs -r -p | wc -l) -ge $MAX_PARALLEL_JOBS ]]; then
        wait -n
      fi
    fi
  done
  
  # Wait for all jobs to complete
  wait
  
  # Generate summary tables after processing all samples
  generate_summary_tables
  
  echo -e "${green}\nPipeline completed successfully!${reset}"
}

# Run main function
main "$@"
