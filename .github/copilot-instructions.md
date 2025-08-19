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
- **Build time: Takes approximately 4-5 seconds. NEVER CANCEL. Set timeout to 60+ seconds.**
- **CRITICAL**: The build will show warnings about deprecated functions (random_shuffle, sprintf) and missing return statements - these are non-critical and do not affect functionality
- Expected build output: Creates `seq2exp` executable (~2.5MB) and intermediate `.o` files (Tools.o, ExprPredictor.o, seq2exp.o)
- Build flags: Uses debug mode by default (`-g -O0`). For release builds, change `CFLAGS = $(DEBUG_FLAGS)` to `CFLAGS = $(RELEASE_FLAGS)` in Makefile
- Requires GSL (GNU Scientific Library) - if build fails with GSL errors, install: `sudo apt-get install libgsl-dev`

### Running the Application
- **CRITICAL**: ALWAYS run from the repository root directory where `iData/` folder is located
- Basic usage with configuration file:
  ```bash
  ./seq2exp -config seq2exp.conf
  ```
- **Runtime: Takes approximately 1-2 seconds for typical datasets with minimal parameters. NEVER CANCEL. Set timeout to 120+ seconds for complex runs with high -nrand values.**
- Command line parameter override example:
  ```bash
  ./seq2exp -config seq2exp.conf -nrand 1 -et 5
  ```
- Traditional command line (without config file):
  ```bash
  ./seq2exp -bs iData/fas.txt -bi iData/topbot2Int.txt -s iData/rhoseq.txt -p iData/synmyout6 -e iData/rhoexp.tab -m iData/factordts.wtmx -f iData/factorexpdts2s.tab -na 1 -i iData/factorinfodts.txt -o BINS -c iData/coopdt.txt -fo out.txt -oo corr -nrand 5 -et 7 -dc iData/coreOPT0dc1.txt -du iData/coreOPT0du1.txt
  ```
- View full usage information: `./seq2exp` (displays all available parameters and model options)

## Validation

### Critical Validation Requirements
**NEVER SKIP**: After ANY code changes, you MUST complete this full validation sequence:
1. **Clean build**: `make clean && make` - verify no fatal errors (warnings are acceptable)
2. **Basic functionality**: `./seq2exp -config seq2exp.conf -nrand 1 -et 5` - verify successful completion
3. **Output verification**: Check `ot.txt` contains "maxBindingWts :" and numerical values
4. **Interface testing**: Test both config file (`-config`) and command line parameter modes
5. **Error handling**: Verify appropriate error messages when files are missing or invalid

### Successful Run Indicators
- Application starts and loads data files without "Cannot open" errors
- Displays matrix data and begins "gradient minimization"
- Shows iteration progress with decreasing objective function values
- Completes with "ending train rng" message
- Exits with status code 0
- Creates output files: `ot.txt`, `ot3.txt`, `pars2.txt`

### Manual Testing Scenario
**ALWAYS run this complete validation sequence after making changes:**
1. **Build**: `make clean && make` (should complete in ~5 seconds without fatal errors)
2. **Execute with minimal parameters**: `./seq2exp -config seq2exp.conf -nrand 1 -et 5` 
3. **Verify successful completion**: Look for "ending train rng" message and exit code 0
4. **Check output files exist**: `ls -la ot*.txt pars*.txt` (should show ot.txt with content, ot3.txt and pars2.txt may be empty)
5. **Validate output content**: `head -5 ot.txt` (should contain "maxBindingWts :" and numerical parameters)
6. **Test traditional CLI**: Run the full command line version to ensure backward compatibility
7. **Test error handling**: Run from wrong directory to verify "Cannot open" error messages

**Expected complete successful run indicators:**
- Application loads data files without "Cannot open" errors
- Displays matrix data (m) and count matrices  
- Shows "Start gradient minimiz minimization" message
- Displays iteration progress with parameter matrices and decreasing objective function values
- Shows final "par _result" and "about to print best par" sections
- Completes with "ending train rng" message
- Exits with status code 0
- Creates `ot.txt` with content like "maxBindingWts :" followed by parameter values

### Troubleshooting
- **Build failures with GSL errors**: Install required package with `sudo apt-get install libgsl-dev`
- **"Cannot open" errors**: Ensure you are running from repository root where `iData/` directory exists
- **Segmentation faults**: Check that all required input files exist in `iData/` directory (15 files should be present)
- **Infinite loops or very long runs**: Use smaller `-nrand` (try 1) and `-et` (try 5) values for testing
- **Empty output files**: Normal for `ot3.txt` and `pars2.txt` with minimal test parameters; `ot.txt` should always contain results
- **Build warnings**: Warnings about deprecated `random_shuffle`, `sprintf`, and missing return statements are non-critical
- **Wrong directory errors**: Must run from `/path/to/MPA/` not from subdirectories or other locations

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
**NEVER SKIP these validation steps after making code modifications:**
1. **Full build and test sequence**:
   ```bash
   make clean && make && ./seq2exp -config seq2exp.conf -nrand 1 -et 5
   ```
   **CRITICAL**: Set timeout to 180+ seconds for this complete sequence. NEVER CANCEL.
2. **Verify output generation**: Check that `ot.txt` contains expected parameter data
3. **Test both interfaces**: Try both config file and traditional command line modes
4. **Performance testing**: For algorithm changes, test with higher `-nrand` values (5-10) to ensure convergence
5. **Error handling**: Test invalid inputs to ensure proper error messages
6. **Cross-platform validation**: Ensure changes work with different compiler versions (warnings are acceptable)

### Performance Tuning
- **Measured timings (Ubuntu 22.04, g++ 13):**
  - Build time: ~4.5 seconds (clean build)
  - Runtime with `-nrand 1 -et 5`: ~1.5 seconds  
  - Runtime with default config (nrand 5): ~3-5 seconds
- Increase `-nrand` (number of random starts) for better convergence but longer runtime (linear scaling)
- Adjust `-et` (energy threshold) to control optimization precision (lower = more precise, longer runtime)
- Use debug build (`-g -O0`) for development, release build (`-O3`) for production performance
- **CRITICAL**: For production runs, expect 5-30 minutes depending on dataset size and parameters. NEVER CANCEL long-running optimization jobs.

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