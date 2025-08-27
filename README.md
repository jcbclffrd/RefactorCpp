# Configuration File Support for seq2exp, built with github copilot coding agent

## Overview

The seq2exp program now supports configuration files in addition to command line arguments. This makes it easier to manage complex parameter sets and run repeated experiments with the same configuration.

## Building from Source

1. Navigate to the source directory: `cd src`
2. Build the executable: `make clean && make`
3. Copy executable to root: `cp seq2exp ../`
4. Return to root directory: `cd ..`
5. Run commands as documented below

## Output Organization

The seq2exp program automatically creates an `oData/` directory for all output files to keep the project organized:

- **Output files**: `ot.txt`, `ot3.txt`, `pars2.txt` - Main results and parameters
- **LaTeX files**: `format.tex` - For generating plots and visualizations  
- **Generated plots**: `plot.pdf` (when LaTeX processing is available)

This separation keeps the project root clean and makes it easy to manage generated data separately from source code and configuration files.

## Usage

### Using a Config File

```bash
./seq2exp -config seq2exp.conf
```

or

```bash
./seq2exp -conf seq2exp.conf
```

### Command Line Override

You can override specific config file parameters with command line arguments:

```bash
./seq2exp -config seq2exp.conf -nrand 10 -et 5
```

### Traditional Command Line (Backward Compatible)

The original command line interface still works:

```bash
./seq2exp -bs iData/fas.txt -bi iData/topbot2Int.txt -s iData/rhoseq.txt -p iData/synmyout6 \
  -e iData/rhoexp.tab -m iData/factordts.wtmx -f iData/factorexpdts2s.tab -na 1 \
  -i iData/factorinfodts.txt -o BINS -c iData/coopdt.txt -fo out.txt -oo corr -nrand 5 \
  -et 7 -dc iData/coreOPT0dc1.txt -du iData/coreOPT0du1.txt
```

## Config File Format

The config file uses a simple key=value format:
- Lines starting with # are comments
- Empty lines are ignored
- Whitespace around keys and values is trimmed

### Example Config File

```
# seq2exp configuration file

# Input files
seqFile = iData/rhoseq.txt
exprFile = iData/rhoexp.tab
motifFile = iData/factordts.wtmx
factorExprFile = iData/factorexpdts2s.tab
factorInfoFile = iData/factorinfodts.txt

# Binding site parameters
bindingSitesFile = iData/fas.txt
bindingIntensityFile = iData/topbot2Int.txt

# Model parameters
paramFile = iData/synmyout6
modelOption = BINS
objOption = corr
coopFile = iData/coopdt.txt

# Training parameters
nAlternations = 1
nRandomStarts = 5
energyThreshold = 7

# Output
outputFile = out.txt
```

## Parameter Mapping

| Config File Key | Command Line | Description |
|----------------|--------------|-------------|
| seqFile | -s | Sequence data file |
| exprFile | -e | Expression data file |
| motifFile | -m | Transcription factor motifs |
| factorExprFile | -f | Factor expression data |
| factorInfoFile | -i | Factor information |
| bindingSitesFile | -bs | Binding site sequences |
| bindingIntensityFile | -bi | Binding intensities |
| paramFile | -p | Parameter file |
| modelOption | -o | Model option |
| objOption | -oo | Objective function |
| coopFile | -c | Cooperativity matrix |
| nAlternations | -na | Number of alternations |
| nRandomStarts | -nrand | Number of random starts |
| energyThreshold | -et | Energy threshold |
| outputFile | -fo | Output file |
| annFile | -a | Annotation file |
| adamiFile | -sa | Adami sequence file |
| exprFile2 | -e2 | Second expression file |
| dcFile | -dc | Dorsal core sequences |
| duFile | -du | Additional motif data |
| repressionFile | -r | Repression file |
| maxContact | -mc | Maximum contact |
| coopDistThr | -ct | Cooperativity distance threshold |
| binwidth | -binwt | Bin width |
| factorIntSigma | -sigma | Sigma parameter |
| repressionDistThr | -rt | Repression distance threshold |
| nExps | -n | Number of experiments |

## Provided Files

- `seq2exp.conf` - Basic configuration file with common parameters