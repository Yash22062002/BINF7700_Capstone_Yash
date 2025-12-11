# SAE Detection Pipeline: 120-Mammal Phylogenetic Analysis
## Computational Pipeline for Synonymous Accelerated Element Detection

**Version:** 2.0  
**Last Updated:** December 2025  
**Author:** Yash Patel (patel.yashm@northeastern.edu)  
**Environment:** Northeastern University Discovery HPC Cluster

---

## Pipeline Overview

This computational pipeline detects Synonymous Accelerated Elements (SAEs) in mammalian protein-coding genes using the FRESCo algorithm. The pipeline processes whole-genome multiple alignment files (MAF format) through gene extraction, preprocessing, evolutionary analysis, and coordinate mapping to produce genome-wide SAE annotations.

### Core Functionality
- Extract gene-level alignments from whole-chromosome MAF files
- Remove reference genome gaps for accurate codon-based analysis
- Detect regions of accelerated synonymous evolution using FRESCo
- Merge overlapping detection windows into discrete genomic elements
- Generate BED format annotations for downstream analysis

---

## System Requirements

### Hardware
- **Memory:** 16GB RAM minimum (32GB recommended for chr1 processing)
- **Storage:** ~500GB for full genome analysis
- **Compute:** SLURM-based HPC cluster with job scheduling

### Software Dependencies
- Python 3.13.5 with pandas and numpy
- R 4.5.2 with tidyverse, ggplot2, knitr, and rmarkdown
- HyPhy 2.5.61 (with FRESCo algorithm)
- MAFfilter 1.3.1
- BEDTools 2.26.0
- dos2unix utility

---

## Pipeline Execution

### Step 1: Data Preparation

Download the 120-mammal alignment files from Hiller Lab repository and prepare gene annotations in GTF format from CCDS database. Extract gene-level alignments using MAFfilter with appropriate coordinate filtering.

### Step 2: Alignment Preprocessing

Remove reference genome gaps using the provided Python script. This step is critical for maintaining proper codon boundaries during evolutionary rate analysis. The script identifies gap positions in the hg38 reference sequence and removes corresponding positions from all 120 species while preserving biological gaps in non-reference species.

### Step 3: SAE Detection with FRESCo

Run FRESCo analysis on individual gene alignments using a 9-codon sliding window approach. The algorithm calculates synonymous substitution rates (dS) and identifies regions with significantly accelerated evolution using Bonferroni-corrected p-values (threshold < 0.05).

### Step 4: Convert to Genomic Coordinates

Transform FRESCo sliding window results into genomic BED coordinates. The conversion script implements codon-overlap based merging to collapse overlapping detection windows into discrete SAE regions. This approach ensures methodological consistency with reference datasets that use non-overlapping regions.

### Step 5: Post-Processing

Clean and standardize output files by removing DOS line endings, sorting coordinates, and merging final results into a single BED file for downstream analysis.

---

## Key Methodological Features

### Window Merging Strategy

The pipeline implements intelligent window merging based on codon overlap rather than simple consecutive window joining. Windows are merged if their codon ranges overlap, creating biologically meaningful SAE regions that accurately represent the underlying evolutionary signal.

### Gap Handling

Reference genome gaps are systematically removed to prevent artifacts in codon-based analysis. This preprocessing step ensures accurate frame preservation and prevents false positive detections caused by alignment artifacts.

### Multi-exon Support

The coordinate mapping system properly handles multi-exon genes, maintaining intronic gaps in BED12 format output while accurately tracking SAE boundaries across splice junctions.

---

## Output Specifications

### FRESCo Output
Tab-delimited text files containing window-level statistics including start/end positions, evolutionary rates (dN/dS), p-values, and SAE/SCE/Normal classifications.

### BED Format
Standard BED12 format with genomic coordinates, gene identifiers (CCDS IDs), strand information, and multi-exon structure. SAE regions are colored red for visualization in genome browsers.

### Mapping Files
Coordinate transformation tables linking alignment positions to genomic coordinates with codon position tracking for accurate frame preservation.

---

## Quality Control

The pipeline includes multiple validation checkpoints:
- Sequence length consistency verification after gap removal
- Overlap detection in final BED output
- Coordinate boundary validation against gene annotations
- Statistical summary generation for detection rates

---

## Performance Considerations

For optimal performance on HPC systems:
- Process chromosomes independently to manage memory usage
- Utilize parallel processing for gene-level analyses
- Implement appropriate SLURM resource allocations
- Clear intermediate files to minimize storage requirements

---

## Troubleshooting Guide

**Common Issues and Solutions:**

**FRESCo alignment errors:** Ensure input alignments have been preprocessed to remove reference gaps

**Memory limitations:** Increase SLURM allocation for large chromosomes or split processing into smaller batches

**Overlapping regions:** Verify usage of v2 scripts with codon-overlap merging functionality

**Module loading failures:** Properly initialize conda environment within SLURM scripts

---

## Citation

When using this pipeline, please cite:

**Software:**
- Hecker & Hiller (2020) - 120-mammal alignment
- Sealfon et al. (2015) - FRESCo algorithm  
- Quinlan & Hall (2010) - BEDTools

**This Pipeline:**
Patel, Y. (2025). SAE Detection Pipeline for 120-Mammal Comparative Genomics. Northeastern University Bioinformatics MS Capstone Project.

---

## Support

For technical questions: patel.yashm@northeastern.edu

For algorithm-specific queries, consult the original publication documentation for FRESCo (Sealfon et al.) and the 120-mammal alignment (Hiller Lab).
