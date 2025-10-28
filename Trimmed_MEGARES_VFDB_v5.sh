#!/bin/bash

# Shell script to automate trimmomatic for multiple samples and subsequent AMR/VF analysis
# Version: 5.0 
# Author: Md. Shaminur Rahman
# Date: August-2025

# Usage
# This script automates running trimmomatic for multiple PE data followed by MEGARES and VFDB analysis
# Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>

# Execution
# Case(1) run on a couple of PE files with extension *.fq
# $ ./Trimmed_MEGARES_VFDB_v5.sh *.fastq.gz

# Case(2) run on a all the fastq files in the current directory (mixed extensions, like .fq. .fastq )
# $ ./Trimmed_MEGARES_VFDB_v5.sh *

# Case(3) skip trimming and run only MEGARES/VFDB analysis on existing paired files
# $ ./Trimmed_MEGARES_VFDB_v5.sh --skip-trim

# Case(4) skip both trimming and MEGARES, run only VFDB analysis on existing paired files
# $ ./Trimmed_MEGARES_VFDB_v5.sh --skip-trim-megares

# Invoke help
# $ ./Trimmed_MEGARES_VFDB_v5.sh -h
# $ ./Trimmed_MEGARES_VFDB_v5.sh --help

##################################################
### CONFIGURATION - EDIT THESE PARAMETERS FIRST ##
##################################################

# System resources
TOTAL_THREADS=80               # Total available CPU threads
THREADS_PER_TRIMMING_JOB=20    # Threads to allocate per trimming job
TOTAL_MEMORY="800G"            # Total memory available for samtools
MEMORY_PER_THREAD="10G"        # Memory per thread for samtools sort

# Calculate maximum memory for samtools sort (threads × memory per thread)
SORT_MEMORY=$(( TOTAL_THREADS * $(echo $MEMORY_PER_THREAD | sed 's/G//') ))G

# Database paths
MEGARES_DB="/home/user/SR/SR_secret/MEGAres/megares_database_v3.00.fasta"
MEGARES_ANNOT="/home/user/SR/SR_secret/MEGAres/megares_to_external_header_mappings_v3.00.csv"
VFDB_DB="/home/user/SR/SR_secret/VFDB/Set-A/Straight-VFDB-Aoutput.fasta"
VFDB_ANNOT="/home/user/SR/SR_secret/VFDB/Set-A/SetA.csv"

# Trimmomatic parameters
ADAPTERS_FILE="/home/user/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa"
HEADCROP=7
LEADING=20
TRAILING=20
SLIDINGWINDOW="20:20"
MINLEN=40

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
Usage: $0 [options] [file_pattern]

Options:
  *.extension       Process files with given extension (fq, fastq, fq.gz, fastq.gz)
  --skip-trim       Skip trimming, use existing files in paired/ directory
  --skip-trim-megares Skip both trimming and MEGARES, run only VFDB analysis
  -h, --help        Show this help message

Configuration (edit script to change):
  Threads: $TOTAL_THREADS total, $THREADS_PER_TRIMMING_JOB per trimming job
  Total memory: $TOTAL_MEMORY
  Memory per thread: $MEMORY_PER_THREAD
  Calculated sort memory: $SORT_MEMORY
  Databases: 
    MEGARES: $MEGARES_DB
    VFDB: $VFDB_DB

Examples:
  $0 *.fq.gz          # Process all gzipped fastq files
  $0 --skip-trim      # Run analysis on pre-trimmed files
${reset}"
  exit 0
}

file_not_found() {
  echo -e "\n${red}FileNotFoundError: One or more input files not found!${reset}"
  echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}\n"
  return 
}

file_name_error() {
  echo -e "\n${red}Filename Error: Paired end file names should contain _R1 _R2${reset}"
  echo -e "${green}Example: test_R1.fq.gz, test_R2.fq.gz${reset}\n"
  return 
}

file_extension_error() {
  echo -e "\n${red}FileExtensionError: Invalid extension for file $1${reset}"
  echo -e "${green}Supported extensions are: <.fq> or <.fastq> or <.fq.gz> or <.fastq.gz>${reset}\n"
  return     
}

create_directories() {
  echo -e "${yellow}Creating required directories...${reset}"
  mkdir -p trimmed_results paired AMR_Results VFDB_Results
}

run_trimmomatic() {
  local R1=$1
  local R2=$2
  local sample_name=$3
  
  R1_pair="paired/${sample_name}_forward_paired.fastq.gz"
  R1_unpair="paired/${sample_name}_forward_unpaired.fastq.gz"
  R2_pair="paired/${sample_name}_reverse_paired.fastq.gz"
  R2_unpair="paired/${sample_name}_reverse_unpaired.fastq.gz"
  
  echo -e "\n${yellow}[Running trimmomatic for sample] ${sample_name}${reset}"
  echo -e "${green}Input files:${reset}"
  echo -e "  R1: $R1"
  echo -e "  R2: $R2"
  echo -e "${green}Output files:${reset}"
  echo -e "  R1 paired: $R1_pair"
  echo -e "  R1 unpaired: $R1_unpair"
  echo -e "  R2 paired: $R2_pair"
  echo -e "  R2 unpaired: $R2_unpair"
  
  trimmomatic PE \
    -threads $THREADS_PER_TRIMMING_JOB \
    -phred33 "$R1" "$R2" "$R1_pair" "$R1_unpair" "$R2_pair" "$R2_unpair" \
    HEADCROP:$HEADCROP \
    ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10:3:TRUE \
    LEADING:$LEADING \
    TRAILING:$TRAILING \
    SLIDINGWINDOW:$SLIDINGWINDOW \
    MINLEN:$MINLEN \
    &>> "trimmed_results/${sample_name}.txt"
  
  # Check if trimming was successful
  if [ $? -ne 0 ]; then
    echo -e "${red}Error: Trimmomatic failed for sample ${sample_name}${reset}"
    return 1
  fi
  
  # Remove unpaired files if they exist
  [ -f "$R1_unpair" ] && rm "$R1_unpair"
  [ -f "$R2_unpair" ] && rm "$R2_unpair"
  
  echo -e "${green}Trimming completed successfully for ${sample_name}${reset}"
}

run_megares() {
  local SAMPLE=$1
  
  echo -e "\n${yellow}[Running MEGARES analysis for sample] ${SAMPLE}${reset}"
  
  # Check if input files exist
  if [ ! -f "paired/${SAMPLE}_forward_paired.fastq.gz" ] || [ ! -f "paired/${SAMPLE}_reverse_paired.fastq.gz" ]; then
    echo -e "${red}Error: Input files for MEGARES analysis not found for sample ${SAMPLE}${reset}"
    return 1
  fi

  # BWA mem with full threads
  echo -e "${green}Running BWA alignment with $TOTAL_THREADS threads...${reset}"
  if ! bwa mem -t $TOTAL_THREADS "$MEGARES_DB" \
    "paired/${SAMPLE}_forward_paired.fastq.gz" \
    "paired/${SAMPLE}_reverse_paired.fastq.gz" \
    > "AMR_Results/${SAMPLE}.amr.alignment.sam"; then
    echo -e "${red}Error: BWA alignment failed for sample ${SAMPLE}${reset}"
    return 1
  fi

  # Convert SAM to BAM
  echo -e "${green}Processing SAM to BAM with $TOTAL_THREADS threads...${reset}"
  if ! samtools view -@ $TOTAL_THREADS -S -b "AMR_Results/${SAMPLE}.amr.alignment.sam" > "AMR_Results/${SAMPLE}.amr.alignment.bam"; then
    echo -e "${red}Error: SAM to BAM conversion failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Sort by queryname for fixmate (required)
  echo -e "${green}Sorting BAM by queryname with $TOTAL_THREADS threads...${reset}"
  if ! samtools sort -@ $TOTAL_THREADS -n -m "$MEMORY_PER_THREAD" "AMR_Results/${SAMPLE}.amr.alignment.bam" -o "AMR_Results/${SAMPLE}.amr.alignment.queryname.bam"; then
    echo -e "${red}Error: BAM sorting by queryname failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Fix mates
  echo -e "${green}Fixing mates with $TOTAL_THREADS threads...${reset}"
  if ! samtools fixmate -@ $TOTAL_THREADS -m "AMR_Results/${SAMPLE}.amr.alignment.queryname.bam" "AMR_Results/${SAMPLE}.amr.alignment.fix.bam"; then
    echo -e "${red}Error: Fixmate failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Sort by coordinate for markdup
  echo -e "${green}Sorting BAM by coordinate with $TOTAL_THREADS threads...${reset}"
  if ! samtools sort -@ $TOTAL_THREADS -m "$MEMORY_PER_THREAD" "AMR_Results/${SAMPLE}.amr.alignment.fix.bam" -o "AMR_Results/${SAMPLE}.amr.alignment.sorted.bam"; then
    echo -e "${red}Error: BAM sorting by coordinate failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Mark duplicates (using markdup instead of rmdup as it's more modern)
  echo -e "${green}Marking duplicates with $TOTAL_THREADS threads...${reset}"
  if ! samtools markdup -@ $TOTAL_THREADS "AMR_Results/${SAMPLE}.amr.alignment.sorted.bam" "AMR_Results/${SAMPLE}.amr.alignment.dedup.bam"; then
    echo -e "${red}Error: Duplicate marking failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Convert back to SAM
  echo -e "${green}Converting back to SAM with $TOTAL_THREADS threads...${reset}"
  if ! samtools view -@ $TOTAL_THREADS -h -o "AMR_Results/${SAMPLE}.amr.alignment.dedup.sam" "AMR_Results/${SAMPLE}.amr.alignment.dedup.bam"; then
    echo -e "${red}Error: BAM to SAM conversion failed for ${SAMPLE}${reset}"
    return 1
  fi

  # Resistome analysis with full threads
  echo -e "${green}Running resistome analysis with $TOTAL_THREADS threads...${reset}"
  if ! resistome \
    -ref_fp "$MEGARES_DB" \
    -annot_fp "$MEGARES_ANNOT" \
    -sam_fp "AMR_Results/${SAMPLE}.amr.alignment.dedup.sam" \
    -gene_fp "AMR_Results/${SAMPLE}.gene.tsv" \
    -group_fp "AMR_Results/${SAMPLE}.group.tsv" \
    -class_fp "AMR_Results/${SAMPLE}.class.tsv" \
    -mech_fp "AMR_Results/${SAMPLE}.mechanism.tsv" \
    -t $TOTAL_THREADS; then
    echo -e "${red}Error: Resistome analysis failed for sample ${SAMPLE}${reset}"
    return 1
  fi
  
  # Clean up intermediate files
  echo -e "${green}Cleaning up intermediate files...${reset}"
  rm -f "AMR_Results/${SAMPLE}.amr.alignment."{sam,bam,queryname.bam,fix.bam,sorted.bam,dedup.sam,dedup.bam}
  
  echo -e "${green}MEGARES analysis completed successfully for ${SAMPLE}${reset}"
}

run_vfdb() {
  local SAMPLE=$1
  
  echo -e "\n${yellow}[Running VFDB analysis for sample] ${SAMPLE}${reset}"
  
  # Check if input files exist
  if [ ! -f "paired/${SAMPLE}_forward_paired.fastq.gz" ] || [ ! -f "paired/${SAMPLE}_reverse_paired.fastq.gz" ]; then
    echo -e "${red}Error: Input files for VFDB analysis not found for sample ${SAMPLE}${reset}"
    return 1
  fi

  # BWA mem with full threads
  echo -e "${green}Running BWA alignment with $TOTAL_THREADS threads...${reset}"
  if ! bwa mem -t $TOTAL_THREADS "$VFDB_DB" \
    "paired/${SAMPLE}_forward_paired.fastq.gz" \
    "paired/${SAMPLE}_reverse_paired.fastq.gz" \
    > "VFDB_Results/${SAMPLE}.vfdb.alignment.sam"; then
    echo -e "${red}Error: BWA alignment failed for sample ${SAMPLE}${reset}"
    return 1
  fi

  # Convert SAM to BAM
  echo -e "${green}Processing SAM to BAM with $TOTAL_THREADS threads...${reset}"
  if ! samtools view -@ $TOTAL_THREADS -S -b "VFDB_Results/${SAMPLE}.vfdb.alignment.sam" > "VFDB_Results/${SAMPLE}.vfdb.alignment.bam"; then
    echo -e "${red}Error: SAM to BAM conversion failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Sort by queryname for fixmate (required)
  echo -e "${green}Sorting BAM by queryname with $TOTAL_THREADS threads...${reset}"
  if ! samtools sort -@ $TOTAL_THREADS -n -m "$MEMORY_PER_THREAD" "VFDB_Results/${SAMPLE}.vfdb.alignment.bam" -o "VFDB_Results/${SAMPLE}.vfdb.alignment.queryname.bam"; then
    echo -e "${red}Error: BAM sorting by queryname failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Fix mates
  echo -e "${green}Fixing mates with $TOTAL_THREADS threads...${reset}"
  if ! samtools fixmate -@ $TOTAL_THREADS -m "VFDB_Results/${SAMPLE}.vfdb.alignment.queryname.bam" "VFDB_Results/${SAMPLE}.vfdb.alignment.fix.bam"; then
    echo -e "${red}Error: Fixmate failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Sort by coordinate for markdup
  echo -e "${green}Sorting BAM by coordinate with $TOTAL_THREADS threads...${reset}"
  if ! samtools sort -@ $TOTAL_THREADS -m "$MEMORY_PER_THREAD" "VFDB_Results/${SAMPLE}.vfdb.alignment.fix.bam" -o "VFDB_Results/${SAMPLE}.vfdb.alignment.sorted.bam"; then
    echo -e "${red}Error: BAM sorting by coordinate failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Mark duplicates (using markdup instead of rmdup as it's more modern)
  echo -e "${green}Marking duplicates with $TOTAL_THREADS threads...${reset}"
  if ! samtools markdup -@ $TOTAL_THREADS "VFDB_Results/${SAMPLE}.vfdb.alignment.sorted.bam" "VFDB_Results/${SAMPLE}.vfdb.alignment.dedup.bam"; then
    echo -e "${red}Error: Duplicate marking failed for ${SAMPLE}${reset}"
    return 1
  fi
  
  # Convert back to SAM
  echo -e "${green}Converting back to SAM with $TOTAL_THREADS threads...${reset}"
  if ! samtools view -@ $TOTAL_THREADS -h -o "VFDB_Results/${SAMPLE}.vfdb.alignment.dedup.sam" "VFDB_Results/${SAMPLE}.vfdb.alignment.dedup.bam"; then
    echo -e "${red}Error: BAM to SAM conversion failed for ${SAMPLE}${reset}"
    return 1
  fi

  # Resistome analysis with full threads
  echo -e "${green}Running resistome analysis with $TOTAL_THREADS threads...${reset}"
  if ! resistome \
    -ref_fp "$VFDB_DB" \
    -annot_fp "$VFDB_ANNOT" \
    -sam_fp "VFDB_Results/${SAMPLE}.vfdb.alignment.dedup.sam" \
    -gene_fp "VFDB_Results/${SAMPLE}.gene.tsv" \
    -group_fp "VFDB_Results/${SAMPLE}.group.tsv" \
    -class_fp "VFDB_Results/${SAMPLE}.class.tsv" \
    -mech_fp "VFDB_Results/${SAMPLE}.mechanism.tsv" \
    -t $TOTAL_THREADS; then
    echo -e "${red}Error: Resistome analysis failed for sample ${SAMPLE}${reset}"
    return 1
  fi
  
  # Clean up intermediate files
  echo -e "${green}Cleaning up intermediate files...${reset}"
  rm -f "VFDB_Results/${SAMPLE}.vfdb.alignment."{sam,bam,queryname.bam,fix.bam,sorted.bam,dedup.sam,dedup.bam}
  
  echo -e "${green}VFDB analysis completed successfully for ${SAMPLE}${reset}"
}

parallel_trim() {
  local files=("$@")
  local pids=()
  local running_jobs=0
  local count=0
  local total_pairs=$(( ${#files[@]} / 2 ))
  local processed_pairs=0

  echo -e "${yellow}Starting parallel trimming with up to $MAX_PARALLEL_JOBS concurrent jobs ($THREADS_PER_TRIMMING_JOB threads each)...${reset}"
  echo -e "${yellow}Total paired samples to process: $total_pairs${reset}"

  for i in "${!files[@]}"; do
    if [[ $((count % 2)) -eq 0 ]]; then
      sample_name=$(basename "${files[$i]}" | awk -F '_R1' '{print $1}')
      R1="${files[$i]}"
      R2="${files[$i]/_R1/_R2}"

      if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo -e "${red}Error: Missing pair for sample ${sample_name}${reset}"
        continue
      fi

      # Wait for a slot if we've reached max parallel jobs
      if [ $running_jobs -ge $MAX_PARALLEL_JOBS ]; then
        wait -n
        running_jobs=$((running_jobs - 1))
      fi

      processed_pairs=$((processed_pairs + 1))
      echo -e "${yellow}Starting trimming for ${sample_name} (job $((running_jobs+1)), pair $processed_pairs/$total_pairs${reset}"
      run_trimmomatic "$R1" "$R2" "$sample_name" &
      pids+=($!)
      running_jobs=$((running_jobs + 1))
    fi
    count=$((count + 1))
  done

  # Wait for all remaining trimming jobs to complete
  echo -e "${yellow}Waiting for remaining $running_jobs trimming jobs to complete...${reset}"
  wait "${pids[@]}"
  
  echo -e "${green}All trimming jobs completed${reset}"
  echo -e "${green}Successfully processed $processed_pairs out of $total_pairs sample pairs${reset}"
}

process_samples() {
  local skip_trim=$1
  local skip_megares=$2
  shift 2
  
  create_directories
  
  if [ "$skip_trim" = false ]; then
    # Verify all files exist before processing
    local missing_files=0
    local valid_files=()
    
    for i in "$@"; do
      if [ ! -f "$i" ]; then
        echo -e "${red}Error: File $i not found${reset}"
        missing_files=$((missing_files + 1))
      elif [[ (${i#*.} == "fastq.gz") || (${i#*.} == "fq.gz") || (${i#*.} == "fastq") || (${i#*.} == "fq") ]]; then
        if echo "$i" | grep -q -e "_R1" -e "_R2"; then
          valid_files+=("$i")
        else
          file_name_error
        fi
      elif [[ ! (${i#*.} == "sh" || ${i#*.} == "sh~") ]]; then
        file_extension_error "$i"
      fi
    done
    
    [ $missing_files -gt 0 ] && file_not_found && return 1
    [ ${#valid_files[@]} -eq 0 ] && echo -e "${red}Error: No valid input files found${reset}" && return 1

    # Run parallel trimming
    parallel_trim "${valid_files[@]}"
  fi
  
  # Run MEGARES and VFDB on all paired samples (sequentially)
  echo -e "\n${yellow}Checking for paired files in paired/ directory...${reset}"
  local found_pairs=0
  local total_pairs=$(ls paired/*_forward_paired.fastq.gz 2>/dev/null | wc -l)
  local current_pair=0
  
  [ $total_pairs -eq 0 ] && echo -e "${yellow}No paired files found in paired/ directory${reset}" && return 0
  
  for paired_file in paired/*_forward_paired.fastq.gz; do
    [ ! -e "$paired_file" ] && continue
    
    local sample_name=$(basename "$paired_file" | awk -F '_forward_paired.fastq.gz' '{print $1}')
    local reverse_file="paired/${sample_name}_reverse_paired.fastq.gz"
    
    if [ -f "$reverse_file" ]; then
      found_pairs=$((found_pairs + 1))
      current_pair=$((current_pair + 1))
      echo -e "\n${green}Processing analysis for sample: ${sample_name} ($current_pair/$total_pairs)${reset}"
      
      # Run MEGARES unless skipped
      if [ "$skip_megares" = false ]; then
        if ! run_megares "$sample_name"; then
          echo -e "${red}MEGARES failed for ${sample_name}, skipping VFDB analysis${reset}"
          continue
        fi
      else
        echo -e "${yellow}Skipping MEGARES analysis for ${sample_name}${reset}"
      fi
      
      # Always run VFDB
      run_vfdb "$sample_name" || echo -e "${red}VFDB failed for ${sample_name}${reset}"
    else
      echo -e "${red}Error: Missing reverse paired file for sample ${sample_name}${reset}"
    fi
  done
  
  if [ $found_pairs -eq 0 ] && [ "$skip_trim" = true ]; then
    echo -e "${red}Error: No valid paired files found in paired/ directory and --skip-trim was specified${reset}"
    return 1
  fi
}

# Main execution
echo -e "\n${yellow}=== Trimmed_MEGARES_VFDB Pipeline ===${reset}"
echo -e "${green}Start time: $(date)${reset}"
echo -e "${yellow}System resources:${reset}"
echo -e "  Total threads: $TOTAL_THREADS"
echo -e "  Threads per trimming job: $THREADS_PER_TRIMMING_JOB"
echo -e "  Max parallel trimming jobs: $MAX_PARALLEL_JOBS"
echo -e "  Total memory: $TOTAL_MEMORY"
echo -e "  Memory per thread: $MEMORY_PER_THREAD"
echo -e "  Calculated sort memory: $SORT_MEMORY\n"

case "$1" in
  -h|--help) usage ;;
  --skip-trim) 
    echo -e "${yellow}Skipping trimming, processing existing paired files${reset}"
    process_samples true false ;;
  --skip-trim-megares)
    echo -e "${yellow}Skipping both trimming and MEGARES, running only VFDB analysis${reset}"
    process_samples true true ;;
  *)
    if [ $# -eq 0 ]; then
      echo -e "${red}Error: No input files specified${reset}"
      usage
      exit 1
    elif [ -e "$1" ]; then
      process_samples false false "$@"
    else
      file_not_found
      exit 1
    fi ;;
esac

# Final status
if [ $? -eq 0 ]; then
  echo -e "\n${green}=== Pipeline completed successfully ===${reset}"
else
  echo -e "\n${red}=== Pipeline completed with errors ===${reset}"
fi

echo -e "${green}End time: $(date)${reset}"
