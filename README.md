# Microbiome Analysis Pipeline
A comprehensive suite of bioinformatics pipelines for metagenomic analysis including taxonomic classification, antimicrobial resistance/virulence factor detection, and functional profiling.

## 🚀 Overview

This repository contains three integrated pipelines for metagenomic data analysis:

1. **Kraken2/Bracken Pipeline** - Taxonomic classification and abundance estimation
2. **AMR/VFDB Analysis Pipeline** - Antimicrobial resistance and virulence factor detection
3. **HUMAnN3 Pipeline** - Functional profiling of microbial communities

## 📋 Requirements

### System Requirements
- **Memory**: Minimum 16GB RAM (recommended 64GB+ for large datasets)
- **Storage**: SSD with sufficient space for databases and temporary files
- **CPU**: Multi-core processor (tested with 80 threads)

### Software Dependencies

#### Core Dependencies (All Pipelines)
- **Bash**: v4.0+
- **Conda**: Miniconda3 or Anaconda
- **Python**: 3.7+
- **Perl**: 5.26+

#### Pipeline-Specific Dependencies

##### Kraken2/Bracken Pipeline
```bash
conda install -c bioconda \
  trimmomatic \
  bbmap \
  kraken2 \
  bracken
```

##### AMR/VFDB Analysis Pipeline
```bash
conda install -c bioconda \
  trimmomatic \
  bwa \
  samtools \
  resistome
```

##### HUMAnN3 Pipeline
```bash
conda install -c biobakery \
  humann \
  metaphlan

conda install -c bioconda \
  trimmomatic \
  bbmap \
  bowtie2 \
  diamond
```

### Database Requirements

#### Kraken2/Bracken
- **Kraken2 Database**: Standard database (~100GB)
- **Host Reference**: Custom host genome (e.g., chicken GRCg6a)

#### AMR/VFDB
- **MEGARES Database**: v3.00 with annotation files
- **VFDB Database**: Set-A with annotation files

#### HUMAnN3
- **UniRef90**: Diamond-formatted database
- **ChocoPhlAn**: pangenome database
- **MetaPhlAn4**: Marker database
- **Utility Mapping**: Database mapping files

## ⚙️ Installation

1. **Clone Repository**
```bash
git clone https://github.com/yourusername/metagenomic-pipelines.git
cd metagenomic-pipelines
```

2. **Set Up Conda Environment**
```bash
# Create and activate environment
conda create -n metagenomics python=3.8
conda activate metagenomics

# Install dependencies (see above for pipeline-specific installations)
```

3. **Configure Database Paths**
Edit the configuration sections in each script to point to your database locations.

## 🛠️ Usage

### 1. Kraken2/Bracken Pipeline

#### Basic Usage
```bash
./Hostfilter_bbmap_Kraken2_Bracken_v6.sh *.fastq.gz
```

#### Options
- `--skip-trim`: Use existing trimmed files
- `-h, --help`: Show help message

#### Example
```bash
# Full analysis
./Hostfilter_bbmap_Kraken2_Bracken_v6.sh sample1_R1.fastq.gz sample1_R2.fastq.gz

# Skip trimming step
./Hostfilter_bbmap_Kraken2_Bracken_v6.sh --skip-trim sample1_R1.fastq.gz
```

### 2. AMR/VFDB Analysis Pipeline

#### Basic Usage
```bash
./Trimmed_MEGARES_VFDB_v5.sh *.fastq.gz
```

#### Options
- `--skip-trim`: Skip trimming, use existing paired files
- `--skip-trim-megares`: Skip trimming and MEGARES, run only VFDB
- `-h, --help`: Show help message

#### Example
```bash
# Full analysis
./Trimmed_MEGARES_VFDB_v5.sh sample1_R1.fastq.gz sample1_R2.fastq.gz

# Skip trimming
./Trimmed_MEGARES_VFDB_v5.sh --skip-trim

# VFDB analysis only
./Trimmed_MEGARES_VFDB_v5.sh --skip-trim-megares
```

### 3. HUMAnN3 Pipeline

#### Basic Usage
```bash
./humann_pipeline.sh <input_dir> <output_dir> <threads_per_sample> <max_parallel_samples>
```

#### Parameters
- `input_dir`: Directory containing paired-end FASTQ files
- `output_dir`: Directory for output files
- `threads_per_sample`: Threads allocated per sample
- `max_parallel_samples`: Maximum number of samples to process in parallel

#### Example
```bash
# Process samples with 16 threads each, max 4 parallel
./humann_pipeline.sh /path/to/input /path/to/output 16 4
```

## 📊 Output Files

### Kraken2/Bracken Pipeline
```
kraken2_results/
├── raw/                    # Raw Kraken2 results
│   └── {sample}/
│       ├── report.tsv
│       └── classifications.txt
├── bracken/
│   ├── species/           # Species-level abundance
│   └── genus/            # Genus-level abundance
├── Final_Species.tsv     # Combined species table
└── Final_Genus.tsv       # Combined genus table
```

### AMR/VFDB Pipeline
```
AMR_Results/
├── {sample}.gene.tsv
├── {sample}.group.tsv
├── {sample}.class.tsv
└── {sample}.mechanism.tsv

VFDB_Results/
├── {sample}.gene.tsv
├── {sample}.group.tsv
├── {sample}.class.tsv
└── {sample}.mechanism.tsv
```

### HUMAnN3 Pipeline
```
output_directory/
├── humann_results/
│   ├── raw_counts/        # Raw gene family counts
│   └── {sample}*          # Processed results
├── trimmed/              # Quality-trimmed files
├── interleaved/          # Interleaved FASTQ files
└── logs/                 # Processing logs
```

## 🔧 Configuration

### Common Parameters (Edit in Scripts)

#### System Resources
```bash
TOTAL_THREADS=80
THREADS_PER_TRIMMING_JOB=20
TOTAL_MEMORY="800G"
```

#### Trimming Parameters
```bash
HEADCROP=7
LEADING=20
TRAILING=20
SLIDINGWINDOW="20:20"
MINLEN=40
```

#### Database Paths
Update these paths in each script's configuration section:
- `KRAKEN_DB`
- `MEGARES_DB` 
- `VFDB_DB`
- `UNIREF90_DB`
- `CHOCOPHLAN_DB`

## 🐛 Troubleshooting

### Common Issues

1. **Memory Errors**
   - Reduce `MAX_PARALLEL_JOBS`
   - Increase `MEMORY_PER_THREAD`
   - Use `--skip-trim` if files are pre-processed

2. **Database Path Errors**
   - Verify all database paths in configuration sections
   - Ensure databases are properly formatted

3. **Missing Dependencies**
   - Run `conda list` to verify installations
   - Check conda channel priorities

### Log Files
Each pipeline generates detailed log files:
- `trimmed_results/` - Trimmomatic logs
- `logs/` - HUMAnN3 processing logs
- Sample-specific logs in output directories

## 📝 Citation

If you use these pipelines in your research, please cite:

Rahman, M. Shaminur, Suborna Islam, Tamanna Jerin Anannya, Adnan Muyeed, Moumita Rahman Sazza, Mohammad Imtiaj Uddin Bhuiyan, and Selina Akter. "Dietary Lees Supplementation Enhances Poultry Health Performance Through Gut Ecosystem Reprogramming." bioRxiv (2026): 2026-01.



## 📞 Support

For questions and support:
- Open an Issue on GitHub
- Contact: s.rahman@just.edu.bd

```

