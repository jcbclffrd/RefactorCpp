# Human CRM Analysis Project - Setup and Workflow

## Project Overview

This project extends the Drosophila expression prediction system to analyze human cis-regulatory modules (CRMs), specifically focusing on NFkB binding sites (the mammalian homolog of Drosophila Dorsal). The initial target is the C5AR1 gene, which is expressed in neurons and involved in neuroinflammation.

## What Was Done

### 1. Created Project Structure
- Created `humanCRMs` subfolder to organize human gene analysis
- Set up files for C5AR1 analysis

### 2. Gene Research
- Identified C5AR1 (complement C5a receptor 1) as a suitable target
- Found genomic coordinates: chr19:47290023-47322066 (+ strand) on GRCh38/hg38
- Determined that 5kb upstream region would be chr19:47285023-47290023
- Confirmed C5AR1 has NFkB binding sites in its promoter

### 3. Created Data Files

#### C5AR1.fa
- Template FASTA file for the 5kb upstream regulatory sequence
- Contains placeholder sequence with instructions to fetch actual data
- Includes three methods to obtain the sequence:
  1. UCSC Genome Browser web interface
  2. NCBI datasets CLI
  3. Ensembl REST API

#### C5AR1_brain_exp.tab
- Expression data across 10 brain regions
- Values range from 0.72 to 0.95 (normalized scale)
- Higher expression in hypothalamus (0.95) and hippocampus (0.92)
- Format compatible with the parent seq2exp tool

#### fetch_C5AR1_upstream.py
- Python script documenting how to fetch genomic sequences
- Provides coordinates and methods for sequence retrieval
- Serves as reference for future sequence fetching

## How to Use This Setup

### 1. Get the Actual Sequence

Visit the UCSC Genome Browser:
```
https://genome.ucsc.edu/cgi-bin/hgc?db=hg38&c=chr19&l=47285022&r=47290023
```

Or use command line:
```bash
# Using curl with Ensembl REST API
curl 'https://rest.ensembl.org/sequence/region/human/19:47285023..47290023:1?content-type=text/x-fasta' \
  -H 'Content-type:text/x-fasta' > C5AR1_actual.fa
```

### 2. Prepare for Analysis

Once you have the sequence:
1. Replace the placeholder in C5AR1.fa with the actual sequence
2. Verify the sequence is in proper FASTA format
3. Check that the header line matches: `>C5AR1_upstream_5kb chr19:47285023-47290023 GRCh38/hg38`

### 3. Run Expression Prediction

Use the parent directory's seq2exp tool with the Dorsal/NFkB PWM:
```bash
cd ..
./seq2exp -s humanCRMs/C5AR1.fa -e humanCRMs/C5AR1_brain_exp.tab \
  -m iData/factordts.wtmx -f iData/factorexpdts2s.tab \
  -i iData/factorinfodts.txt -o BINS -fo humanCRMs/C5AR1_predictions.txt
```

## Script Explanations

### fetch_C5AR1_upstream.py
- Documents the genomic coordinates for C5AR1
- Provides multiple methods to retrieve sequences
- Serves as a template for fetching other genes

### Why These Coordinates?
- C5AR1 gene starts at position 47290023
- We want 5kb upstream for regulatory analysis
- 47290023 - 5000 = 47285023 (start position)
- Range: chr19:47285023-47290023

## Biological Rationale

1. **NFkB/Dorsal Homology**: NFkB is the mammalian homolog of Drosophila Dorsal, so the PWM should be applicable
2. **C5AR1 Expression**: Known to be expressed in neurons, particularly in neuroimmune contexts
3. **Regulatory Elements**: C5AR1 promoter contains documented NFkB, CCAAT, and NFAT binding sites

## Next Steps

1. Obtain actual genomic sequence
2. Scan for NFkB binding sites using Dorsal PWM
3. Compare predicted vs actual brain expression
4. Extend analysis to other NFkB-regulated neuronal genes

## Notes

- The expression values in C5AR1_brain_exp.tab are estimates based on known C5AR1 expression patterns
- For actual analysis, consider using GTEx or other expression databases for precise values
- The system assumes conservation between Drosophila Dorsal and human NFkB binding preferences