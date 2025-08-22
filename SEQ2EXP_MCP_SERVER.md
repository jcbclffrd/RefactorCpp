# seq2exp MCP Server Documentation

## Overview

The seq2exp MCP Server provides a clean API interface to the real seq2exp bioinformatics application for gene expression prediction. It accepts configuration file content directly as text and executes the actual seq2exp C++ application to generate predictions and PDF plots.

## Key Features

- **Direct Config File Interface**: Accepts seq2exp configuration as plain text (no complex JSON transformation)
- **Real Application Execution**: Runs the actual seq2exp C++ bioinformatics application
- **Fast Execution**: ~1.5 seconds with minimal test parameters
- **PDF Generation**: Automatically processes LaTeX output to generate plot.pdf
- **Comprehensive Error Handling**: Proper validation, timeouts, and error reporting
- **RESTful API**: Standard HTTP endpoints with OpenAPI documentation

## Quick Start

### 1. Start the Server

```bash
cd /home/runner/work/MPA/MPA
python3 seq2exp_mcp_server.py
```

Server will be available at `http://localhost:8083`

### 2. Get Configuration Template

```bash
curl -X GET "http://localhost:8083/tools/seq2exp_config_template"
```

### 3. Run Prediction

```bash
curl -X POST "http://localhost:8083/tools/seq2exp_predict" \
  -H "Content-Type: application/json" \
  -d '{
    "config_content": "# Fast test config\nseqFile = iData/rhoseq.txt\nexprFile = iData/rhoexp.tab\nmotifFile = iData/factordts.wtmx\nfactorExprFile = iData/factorexpdts2s.tab\nfactorInfoFile = iData/factorinfodts.txt\nbindingSitesFile = iData/fas.txt\nbindingIntensityFile = iData/topbot2Int.txt\nparamFile = iData/synmyout6\nmodelOption = BINS\nobjOption = corr\ncoopFile = iData/coopdt.txt\nnAlternations = 1\nnRandomStarts = 1\nenergyThreshold = 5\noutputFile = out.txt\ndcFile = iData/coreOPT0dc1.txt\nduFile = iData/coreOPT0du1.txt"
  }'
```

## API Endpoints

### 1. Configuration Template

**GET** `/tools/seq2exp_config_template`

Returns the default seq2exp configuration file as a template.

**Response:**
```json
{
  "template": "# Complete seq2exp configuration file\n...",
  "description": "Default seq2exp configuration template...",
  "usage": "Copy this template and modify..."
}
```

### 2. Configuration Validation

**POST** `/tools/seq2exp_validate_config`

Validates configuration file content without running seq2exp.

**Request Body:**
```json
{
  "config_content": "seqFile = iData/rhoseq.txt\nexprFile = iData/rhoexp.tab\n..."
}
```

**Response:**
```json
{
  "valid": true,
  "issues": [],
  "parsed_config": {
    "seqFile": "iData/rhoseq.txt",
    "exprFile": "iData/rhoexp.tab",
    ...
  },
  "parameter_count": 16
}
```

### 3. Gene Expression Prediction

**POST** `/tools/seq2exp_predict`

Main prediction endpoint that executes seq2exp and returns results.

**Request Body:**
```json
{
  "config_content": "# seq2exp configuration\nseqFile = iData/rhoseq.txt\n..."
}
```

**Response:**
```json
{
  "success": true,
  "execution_time": 1.544,
  "results": {
    "main_output": " maxBindingWts : \n0.891\t1.000\t0.990\t\n...",
    "binding_weights": [0.891, 1.0, 0.99],
    "objective_function": 14.962
  },
  "pdf_path": null,
  "stdout": "m \n0.000000 1.000000 62.000000...",
  "config_content": "# Fast test config\nseqFile = iData/rhoseq.txt\n...",
  "returncode": 0
}
```

### 4. Get Latest Results

**GET** `/tools/seq2exp_get_results`

Retrieves results from the latest seq2exp run without re-executing.

**Response:**
```json
{
  "results": {
    "main_output": " maxBindingWts : \n0.891\t1.000\t0.990\t\n...",
    "binding_weights": [0.891, 1.0, 0.99],
    "objective_function": 14.962
  },
  "pdf_available": true,
  "pdf_url": "/results/pdf/plot.pdf",
  "timestamp": 1692123456.789
}
```

### 5. PDF Download

**GET** `/results/pdf/{filename}`

Downloads generated PDF files (e.g., plot.pdf).

## Configuration Parameters

The server accepts the same configuration parameters as the original seq2exp application. Key parameters:

### Required Parameters
- `seqFile`: Input sequence data (FASTA format)
- `exprFile`: Expression data (tab-delimited)
- `motifFile`: Transcription factor motifs
- `factorExprFile`: Factor expression data
- `outputFile`: Output file for predictions

### Performance Parameters
- `nRandomStarts`: Number of random starts (1 for fast testing, 5+ for accuracy)
- `energyThreshold`: Energy threshold (5 for fast testing, 7+ for accuracy)
- `nAlternations`: Number of alternations (1 for training)

### Model Options
- `modelOption`: BINS, Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited
- `objOption`: corr, SSE, Cross_Corr

## Usage Examples

### Example 1: Fast Testing Configuration

```bash
curl -X POST "http://localhost:8083/tools/seq2exp_predict" \
  -H "Content-Type: application/json" \
  -d '{
    "config_content": "# Fast test configuration\nseqFile = iData/rhoseq.txt\nexprFile = iData/rhoexp.tab\nmotifFile = iData/factordts.wtmx\nfactorExprFile = iData/factorexpdts2s.tab\nfactorInfoFile = iData/factorinfodts.txt\nbindingSitesFile = iData/fas.txt\nbindingIntensityFile = iData/topbot2Int.txt\nparamFile = iData/synmyout6\nmodelOption = BINS\nobjOption = corr\ncoopFile = iData/coopdt.txt\nnAlternations = 1\nnRandomStarts = 1\nenergyThreshold = 5\noutputFile = out.txt\ndcFile = iData/coreOPT0dc1.txt\nduFile = iData/coreOPT0du1.txt"
  }'
```

**Expected Response Time:** ~1.5 seconds

### Example 2: Production Configuration

```bash
curl -X POST "http://localhost:8083/tools/seq2exp_predict" \
  -H "Content-Type: application/json" \
  -d '{
    "config_content": "# Production configuration\nseqFile = iData/rhoseq.txt\nexprFile = iData/rhoexp.tab\nmotifFile = iData/factordts.wtmx\nfactorExprFile = iData/factorexpdts2s.tab\nfactorInfoFile = iData/factorinfodts.txt\nbindingSitesFile = iData/fas.txt\nbindingIntensityFile = iData/topbot2Int.txt\nparamFile = iData/synmyout6\nmodelOption = BINS\nobjOption = corr\ncoopFile = iData/coopdt.txt\nnAlternations = 1\nnRandomStarts = 5\nenergyThreshold = 7\noutputFile = out.txt\ndcFile = iData/coreOPT0dc1.txt\nduFile = iData/coreOPT0du1.txt"
  }'
```

**Expected Response Time:** ~3-5 seconds

### Example 3: Configuration Validation

```bash
curl -X POST "http://localhost:8083/tools/seq2exp_validate_config" \
  -H "Content-Type: application/json" \
  -d '{
    "config_content": "seqFile = iData/rhoseq.txt\nexprFile = iData/rhoexp.tab\ninvalidParameter = badValue"
  }'
```

## Python Client Example

```python
import aiohttp
import asyncio
import json

async def run_seq2exp_prediction():
    config = """# Fast test config
seqFile = iData/rhoseq.txt
exprFile = iData/rhoexp.tab
motifFile = iData/factordts.wtmx
factorExprFile = iData/factorexpdts2s.tab
factorInfoFile = iData/factorinfodts.txt
bindingSitesFile = iData/fas.txt
bindingIntensityFile = iData/topbot2Int.txt
paramFile = iData/synmyout6
modelOption = BINS
objOption = corr
coopFile = iData/coopdt.txt
nAlternations = 1
nRandomStarts = 1
energyThreshold = 5
outputFile = out.txt
dcFile = iData/coreOPT0dc1.txt
duFile = iData/coreOPT0du1.txt"""

    async with aiohttp.ClientSession() as session:
        async with session.post(
            'http://localhost:8083/tools/seq2exp_predict',
            json={'config_content': config}
        ) as response:
            result = await response.json()
            
            if result['success']:
                print(f"Prediction completed in {result['execution_time']:.2f}s")
                print(f"Binding weights: {result['results']['binding_weights']}")
                print(f"Objective function: {result['results']['objective_function']}")
            else:
                print(f"Prediction failed: {result['error']}")

# Run the example
asyncio.run(run_seq2exp_prediction())
```

## Error Handling

The server provides comprehensive error handling:

### Common Error Scenarios

1. **Missing Input Files**
   ```json
   {
     "valid": false,
     "issues": ["File not found: /path/to/missing/file.txt"]
   }
   ```

2. **seq2exp Execution Failure**
   ```json
   {
     "success": false,
     "error": "seq2exp failed: Error message from application",
     "returncode": 1
   }
   ```

3. **Timeout**
   ```json
   {
     "success": false,
     "error": "seq2exp execution timed out after 120 seconds"
   }
   ```

4. **Invalid Configuration Format**
   ```json
   {
     "valid": false,
     "issues": ["Line 5: Missing '=' - invalid line format"]
   }
   ```

## Performance Notes

- **Testing Parameters**: Use `nRandomStarts=1, energyThreshold=5` for ~1-2 second execution
- **Production Parameters**: Use `nRandomStarts=5+, energyThreshold=7+` for accurate results (3-30 minutes)
- **Concurrent Requests**: Server handles multiple simultaneous requests
- **Timeout**: Default 120-second timeout for seq2exp execution

## OpenAPI Documentation

When the server is running, visit `http://localhost:8083/docs` for interactive API documentation with Swagger UI.

## Dependencies

### System Requirements
- Ubuntu/Debian Linux
- GSL (GNU Scientific Library): `sudo apt-get install libgsl-dev`
- LaTeX (for PDF generation): `sudo apt-get install texlive`

### Python Requirements
- FastAPI
- uvicorn
- aiofiles
- python-multipart

Install with:
```bash
pip install fastapi uvicorn aiofiles python-multipart
```

## Architecture

```
Config Content (Text) → Temp Config File → ./seq2exp → Results + format.tex → scriptse.sh → plot.pdf
```

The server follows the problem statement requirement of accepting config file content directly as text, creating temporary config files, and executing the real seq2exp application without reimplementing its functionality.