# Copilot Instructions for MPA Codebase

## Project Overview
This repository centers on the `seq2exp` program, which predicts gene expression from sequence and motif data. The codebase is primarily C++ and uses configuration files for flexible experiment setup.

## Architecture & Key Components
- **Main Programs:**
  - `seq2exp.cpp`: Core logic for running experiments and parsing config/command-line arguments.
  - `ExprPredictor.cpp/h`: Implements the main prediction model and related algorithms.
  - `Tools.cpp/h`: Utility functions for file I/O, parsing, and data manipulation.
- **Data Directory:**
  - `iData/`: Contains all input data files (sequences, motifs, expression, etc.) used in experiments.
- **Configuration:**
  - `seq2exp.conf`: Key-value config file for experiment parameters. See README for mapping to CLI args.

## Developer Workflows
- **Build:**
  - Use the provided `Makefile` to build all binaries: `make`
- **Run Experiments:**
  - Typical usage: `./seq2exp -config seq2exp.conf`
  - Command-line arguments override config file values.
- **Debugging:**
  - Add debug output via `std::cerr` or use gdb/lldb for step-through debugging.
- **Data Files:**
  - All input/output files referenced in configs/CLI are expected in `iData/` or workspace root.

## Project-Specific Conventions
- **Config Parsing:**
  - Config files use `key = value` format, with comments starting with `#`.
  - Whitespace is trimmed; empty lines ignored.
- **Parameter Mapping:**
  - See README for mapping between config keys and CLI flags.
- **Backward Compatibility:**
  - Direct CLI usage is supported for legacy workflows.
- **File Naming:**
  - Input files are named for their biological role (e.g., `fas.txt`, `rhoexp.tab`).

## Integration & Dependencies
- No external libraries required; pure C++.
- All dependencies are local files and standard C++ STL.

## Examples
- **Run with config:**
  ```bash
  ./seq2exp -config seq2exp.conf
  ```
- **Override config values:**
  ```bash
  ./seq2exp -config seq2exp.conf -nrand 10 -et 5
  ```
- **Build:**
  ```bash
  make
  ```

## Key Files & Directories
- `seq2exp.cpp`, `ExprPredictor.cpp/h`, `Tools.cpp/h`: Main source files
- `iData/`: Input data files
- `seq2exp.conf`: Example config
- `Makefile`: Build instructions
- `README.md`: Parameter mapping and usage details

---

If any conventions or workflows are unclear, please provide feedback so this guide can be improved.
