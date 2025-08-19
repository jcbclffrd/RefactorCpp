# seq2exp Bioinformatics Application

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

seq2exp is a C++ bioinformatics application for sequence-to-expression prediction using machine learning optimization techniques. It processes DNA sequences, expression data, and transcription factor motifs to predict gene expression patterns.

## Working Effectively

### Dependencies and Setup
- Install required system packages:
  ```bash
  sudo apt-get update && sudo apt-get install -y libgsl-dev build-essential
  ```
- The application requires GSL (GNU Scientific Library) for mathematical computations
- Uses standard g++ compiler with C++13 features

### Build Process
- Clean and build the application:
  ```bash
  make clean && make
  ```
- Build time: Takes approximately 5 seconds. NEVER CANCEL. Set timeout to 30+ seconds.
- The build produces warnings about deprecated functions but these are non-critical and do not affect functionality
- Expected build output: Creates `seq2exp` executable and intermediate `.o` files
- Build flags: Uses debug mode by default (`-g -O0`). For release builds, change to `RELEASE_FLAGS` in Makefile

### Running the Application
- ALWAYS run from the repository root directory where `iData/` folder is located
- Basic usage with configuration file:
  ```bash
  ./seq2exp -config seq2exp.conf
  ```
- Runtime: Takes approximately 1-2 seconds for typical datasets. NEVER CANCEL. Set timeout to 60+ seconds for large datasets.
- Command line parameter override example:
  ```bash
  ./seq2exp -config seq2exp.conf -nrand 1 -et 5
  ```
- Traditional command line (without config file):
  ```bash
  ./seq2exp -bs iData/fas.txt -bi iData/topbot2Int.txt -s iData/rhoseq.txt -p iData/synmyout6 -e iData/rhoexp.tab -m iData/factordts.wtmx -f iData/factorexpdts2s.tab -na 1 -i iData/factorinfodts.txt -o BINS -c iData/coopdt.txt -fo out.txt -oo corr -nrand 5 -et 7 -dc iData/coreOPT0dc1.txt -du iData/coreOPT0du1.txt
  ```

## Validation

### Successful Run Indicators
- Application starts and loads data files without "Cannot open" errors
- Displays matrix data and begins "gradient minimization"
- Shows iteration progress with decreasing objective function values
- Completes with "ending train rng" message
- Exits with status code 0
- Creates output files: `ot.txt`, `ot3.txt`, `pars2.txt`

### Manual Testing Scenario
After making changes, ALWAYS test by running:
1. Build: `make clean && make`
2. Execute: `./seq2exp -config seq2exp.conf -nrand 1 -et 5`
3. Verify output files exist: `ls -la ot*.txt pars*.txt`
4. Check output file content: `head -5 ot.txt`
5. Expected content includes lines like "maxBindingWts :" and numerical parameters

### Troubleshooting
- If build fails with GSL errors: Install `libgsl-dev` package
- If "Cannot open" errors: Ensure running from repository root where `iData/` exists
- If segmentation faults: Check that all required input files exist in `iData/`
- If infinite loops: Use smaller `-nrand` and `-et` values for testing

## Configuration

### Config File Format
- Uses `seq2exp.conf` for default parameters
- Key=value format with # for comments
- Critical parameters:
  - `seqFile`, `exprFile`, `motifFile`, `factorExprFile`, `outputFile` - required input/output files
  - `nRandomStarts`, `energyThreshold` - control optimization
  - `modelOption` - algorithm choice (BINS, Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited)

### Command Line Options
- Use `-config` or `-conf` to specify config file
- Any config parameter can be overridden via command line using its flag
- Use `./seq2exp` without arguments to see full usage information

## Key Projects and Files

### Source Code Structure
- `seq2exp.cpp` - Main application entry point, command line parsing, configuration loading
- `ExprPredictor.cpp/.h` - Core prediction algorithms and machine learning optimization  
- `Tools.cpp/.h` - Mathematical utilities, matrix operations, GSL integration
- `Makefile` - Build configuration with GSL library linking

### Data Files (iData/ directory)
- `rhoseq.txt` - Input sequence data (FASTA format)
- `rhoexp.tab` - Expression data (tab-delimited)
- `factordts.wtmx` - Transcription factor motifs
- `factorexpdts2s.tab` - Factor expression data
- `factorinfodts.txt` - Factor information (activator/repressor roles)
- `fas.txt` - Binding site sequences
- `topbot2Int.txt` - Binding intensities
- Additional core sequence files: `coreOPT0dc1.txt`, `coreOPT0du1.txt`

### Configuration and Documentation
- `seq2exp.conf` - Complete example configuration file with all parameters
- `README.md` - Usage documentation and parameter mapping
- `scriptse.sh` - LaTeX processing script for generating plots

## Common Tasks

### Testing Changes
1. Always build and test after modifications:
   ```bash
   make clean && make && ./seq2exp -config seq2exp.conf -nrand 1 -et 5
   ```
2. For debugging, run with minimal parameters to reduce execution time
3. Check that output files are generated and contain expected numerical results

### Performance Tuning
- Increase `-nrand` (number of random starts) for better convergence but longer runtime
- Adjust `-et` (energy threshold) to control optimization precision
- Use debug build (`-g -O0`) for development, release build (`-O3`) for production

### Repository Root Contents
```
.
├── .github/               # GitHub configuration
├── iData/                # Required input data files
├── ExprPredictor.cpp     # Core prediction algorithms 
├── ExprPredictor.h       # Algorithm headers
├── Tools.cpp             # Mathematical utilities
├── Tools.h               # Utility headers  
├── seq2exp.cpp           # Main application
├── Makefile              # Build configuration
├── seq2exp.conf          # Example configuration
├── README.md             # Documentation
├── scriptse.sh           # LaTeX processing
├── plot.tex/.pdf         # Generated plots
└── format.tex            # LaTeX formatting
```

### Build Artifacts (Not Committed)
- `seq2exp` - Main executable
- `*.o` - Object files  
- `ot*.txt`, `pars*.txt` - Output files from runs