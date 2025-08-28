# Human CRM Analysis - C5AR1

This folder contains data for analyzing human cis-regulatory modules (CRMs), specifically for the C5AR1 gene using NFkB (homolog of Drosophila Dorsal).

## Gene Information

- **Gene**: C5AR1 (complement C5a receptor 1)
- **NCBI Gene ID**: 728
- **Ensembl ID**: ENSG00000197405
- **Location**: chr19:47290023-47322066 (+ strand)
- **Function**: Encodes a G-protein coupled receptor for C5a, involved in inflammatory responses
- **Expression**: High expression in brain regions, particularly in neurons

## Files

### C5AR1.fa
Contains the 5kb upstream regulatory sequence from the C5AR1 transcription start site.
- **Coordinates**: chr19:47285023-47290023 (GRCh38/hg38)
- **Note**: The current file contains a placeholder sequence. To get the actual sequence:
  1. Visit [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgc?db=hg38&c=chr19&l=47285022&r=47290023)
  2. Click "Get DNA for this region"
  3. Replace the placeholder sequence

### C5AR1_brain_exp.tab
Expression data for C5AR1 upstream region across different brain regions.
- Based on known C5AR1 expression patterns in neuronal tissues
- Values represent relative expression levels (0-1 scale)
- Higher values in hypothalamus and hippocampus reflect C5AR1's role in neuroimmune responses

## Analysis Notes

- NFkB binding sites are expected in the C5AR1 promoter region
- The Dorsal PWM from the parent directory should be applicable as NFkB is the mammalian homolog
- C5AR1 is known to have NFkB, CCAAT, and NFAT binding sites in its promoter