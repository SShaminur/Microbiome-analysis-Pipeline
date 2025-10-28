# Microbiome Analysis Pipeline

A comprehensive suite of bash and R scripts for metagenomic analysis including taxonomic profiling, AMR detection, and functional analysis.

## Features

- **Quality Control**: Trimmomatic and BBduk preprocessing
- **Taxonomic Profiling**: Kraken2/Bracken with host filtering
- **Antimicrobial Resistance**: MEGARES database analysis
- **Virulence Factors**: VFDB database screening
- **Functional Profiling**: HUMAnN3 pathway analysis


## Quick Start

```bash
# Clone repository
git clone https://github.com/yourusername/microbiome-pipeline.git
cd microbiome-pipeline

# Run taxonomic profiling
./scripts/taxonomic-profiling/Hostfilter_bbmap_Kraken2_Bracken_v6.sh *.fastq.gz

# Run AMR/VF analysis
./scripts/amr-analysis/Trimmed_MEGARES_VFDB_v5.sh *.fastq.gz
