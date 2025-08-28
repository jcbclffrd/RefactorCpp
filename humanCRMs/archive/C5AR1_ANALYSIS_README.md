# C5AR1 Expression Analysis Setup

## Overview

This setup allows running seq2exp analysis on the human C5AR1 gene using brain region expression data and the transcription factors dl (dorsal), tw (twist), and sn (snail).

## Directory Structure

```
humanCRMs/
├── C5AR1.fa                       # C5AR1 upstream sequence (5kb)
├── C5AR1_brain_exp.tab           # Expression data across 10 brain regions
├── c5ar1_factorexp_brain.tab    # TF expression levels in brain regions
├── c5ar1_data/                   # Supporting data files
│   ├── factordts.wtmx            # TF binding motifs (dl, tw, sn)
│   ├── factorinfodts.txt         # TF information
│   ├── coopdt.txt                # Cooperativity matrix
│   ├── synmyout6                 # Initial parameters
│   ├── coreOPT0dc1.txt          # Dorsal core sequences
│   ├── coreOPT0du1.txt          # Additional motif data
│   ├── c5ar1_binding_sites.txt  # Predicted binding sites
│   └── c5ar1_binding_intensity.txt # Binding intensities
└── c5ar1_predictions.txt         # Output predictions

c5ar1_seq2exp.conf                # Configuration file
run_c5ar1_analysis.sh             # Run script
```

## Brain Regions Analyzed

The expression data includes 10 brain regions:
1. Frontal Cortex
2. Hippocampus
3. Cerebellum
4. Striatum
5. Hypothalamus
6. Amygdala
7. Thalamus
8. Brain Stem
9. Pituitary
10. Spinal Cord

## Transcription Factors

The analysis uses three key transcription factors:
- **dl (dorsal)**: Nuclear factor involved in dorsal-ventral patterning
- **tw (twist)**: Basic helix-loop-helix transcription factor
- **sn (snail)**: Zinc finger transcription factor

## Running the Analysis

1. Make sure seq2exp is compiled:
   ```bash
   make seq2exp
   ```

2. Run the analysis:
   ```bash
   ./run_c5ar1_analysis.sh
   ```

   Or directly:
   ```bash
   ./seq2exp -config c5ar1_seq2exp.conf
   ```

## Configuration

The `c5ar1_seq2exp.conf` file contains all parameters for the analysis:
- Input sequence: C5AR1 5kb upstream region
- Expression data: Brain region-specific expression levels
- TF motifs and expression: dl, tw, sn factors
- Model parameters: BINS model with correlation objective

## Output

Results are saved to `humanCRMs/c5ar1_predictions.txt` containing:
- Predicted expression values for each brain region
- Model performance metrics
- Trained parameter values

## Notes

- The binding sites and intensities are placeholder values for initial analysis
- TF expression levels in brain regions are estimated based on typical neural expression patterns
- The original motif and cooperativity data from the Drosophila model are retained