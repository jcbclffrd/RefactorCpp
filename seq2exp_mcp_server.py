#!/usr/bin/env python3
"""
seq2exp MCP Server - Real bioinformatics application wrapper
This server provides MCP tools that directly execute the seq2exp application
and process its LaTeX output to generate PDF plots.
"""

import asyncio
import os
import tempfile
import time
from pathlib import Path
from typing import Dict, Any, Optional

import aiofiles
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field


# Configuration
WORKING_DIR = "/home/runner/work/MPA/MPA"
SEQ2EXP_EXECUTABLE = f"{WORKING_DIR}/seq2exp"
SCRIPT_SE_PATH = f"{WORKING_DIR}/scriptse.sh"
DEFAULT_CONFIG_PATH = f"{WORKING_DIR}/seq2exp.conf"


class Seq2ExpRequest(BaseModel):
    config_content: str = Field(
        description="Complete seq2exp configuration file content",
        example="""# seq2exp configuration
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
    )


class Seq2ExpConfigManager:
    """Handles seq2exp configuration file validation and management"""
    
    @staticmethod
    def validate_config_content(config_content: str) -> Dict[str, Any]:
        """Validate config file content and return status"""
        issues = []
        config_dict = {}
        
        # Parse config content
        for line_num, line in enumerate(config_content.split('\n'), 1):
            line = line.strip()
            if line and not line.startswith('#'):
                if '=' in line:
                    try:
                        key, value = line.split('=', 1)
                        config_dict[key.strip()] = value.strip()
                    except ValueError:
                        issues.append(f"Line {line_num}: Invalid format - {line}")
                else:
                    issues.append(f"Line {line_num}: Missing '=' - {line}")
        
        # Check for required files in iData directory
        file_keys = [k for k in config_dict.keys() if k.endswith('File')]
        for key in file_keys:
            if key in config_dict:
                filepath = f"{WORKING_DIR}/{config_dict[key]}"
                if not os.path.exists(filepath):
                    issues.append(f"File not found: {filepath} (for {key})")
        
        # Check for required parameters
        required_params = ['seqFile', 'exprFile', 'motifFile', 'factorExprFile', 'outputFile']
        for param in required_params:
            if param not in config_dict:
                issues.append(f"Missing required parameter: {param}")
        
        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "parsed_config": config_dict,
            "parameter_count": len(config_dict)
        }
    
    @staticmethod
    async def get_default_config() -> str:
        """Return default seq2exp.conf content as template"""
        try:
            async with aiofiles.open(DEFAULT_CONFIG_PATH, 'r') as f:
                return await f.read()
        except FileNotFoundError:
            # Fallback template if default config doesn't exist
            return """# seq2exp configuration file
# This file contains all available parameters
# Lines starting with # are comments
# Uncomment and modify parameters as needed

# === REQUIRED INPUT FILES ===
seqFile = iData/rhoseq.txt                    # -s: Sequence data (FASTA format)
exprFile = iData/rhoexp.tab                   # -e: Expression data (tab-delimited)
motifFile = iData/factordts.wtmx              # -m: Transcription factor motifs
factorExprFile = iData/factorexpdts2s.tab     # -f: Factor expression data
outputFile = out.txt                          # -fo: Output file for predictions

# === BINDING SITE PARAMETERS ===
bindingSitesFile = iData/fas.txt              # -bs: Binding site sequences (FASTA format)
bindingIntensityFile = iData/topbot2Int.txt   # -bi: Binding intensities (tab-delimited)

# === MODEL PARAMETERS ===
paramFile = iData/synmyout6                   # -p: Parameter file for initial values
modelOption = BINS                            # -o: Model option (BINS, Logistic, Direct, Quenching, ChrMod_Unlimited, ChrMod_Limited)
objOption = corr                              # -oo: Objective function (corr, SSE, Cross_Corr)
coopFile = iData/coopdt.txt                   # -c: Cooperativity matrix

# === REQUIRED FOR CERTAIN MODELS ===
factorInfoFile = iData/factorinfodts.txt      # -i: Factor information (activator/repressor roles)

# === TRAINING PARAMETERS ===
nAlternations = 1                             # -na: Number of alternations (1 for training, 0 for testing)
nRandomStarts = 5                             # -nrand: Number of random starts
energyThreshold = 7                           # -et: Energy threshold

# === ADDITIONAL DATA FILES ===
dcFile = iData/coreOPT0dc1.txt                # -dc: Dorsal core sequences
duFile = iData/coreOPT0du1.txt                # -du: Additional motif data"""


class Seq2ExpExecutor:
    """Handles seq2exp execution and output processing"""
    
    def __init__(self, working_dir: str = WORKING_DIR):
        self.working_dir = working_dir
        self.executable_path = SEQ2EXP_EXECUTABLE
        
    async def run_prediction(self, config_content: str) -> Dict[str, Any]:
        """Main prediction execution"""
        # 1. Create temporary config file from content
        with tempfile.NamedTemporaryFile(mode='w', suffix='.conf', delete=False, dir=self.working_dir) as f:
            f.write(config_content)
            temp_config_path = f.name
            
        try:
            # 2. Execute seq2exp
            start_time = time.time()
            process = await asyncio.create_subprocess_exec(
                './seq2exp', '-config', os.path.basename(temp_config_path),
                cwd=self.working_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await asyncio.wait_for(
                process.communicate(), timeout=120  # 2 minute timeout
            )
            
            execution_time = time.time() - start_time
            
            if process.returncode != 0:
                return {
                    "success": False,
                    "error": f"seq2exp failed: {stderr.decode()}",
                    "returncode": process.returncode,
                    "stdout": stdout.decode(),
                    "execution_time": execution_time
                }
            
            # 3. Parse output files
            results = await self._parse_seq2exp_output()
            
            # 4. Generate PDF if format.tex exists
            pdf_path = await self._generate_pdf()
            
            return {
                "success": True,
                "execution_time": execution_time,
                "results": results,
                "pdf_path": pdf_path,
                "stdout": stdout.decode(),
                "config_content": config_content,
                "returncode": process.returncode
            }
            
        except asyncio.TimeoutError:
            return {
                "success": False,
                "error": "seq2exp execution timed out after 120 seconds",
                "execution_time": 120.0
            }
        except Exception as e:
            return {
                "success": False,
                "error": f"Execution error: {str(e)}",
                "execution_time": time.time() - start_time if 'start_time' in locals() else 0
            }
        finally:
            # Cleanup temp config file
            try:
                os.unlink(temp_config_path)
            except:
                pass
    
    async def _parse_seq2exp_output(self) -> Dict[str, Any]:
        """Parse seq2exp output files"""
        results = {}
        
        # Parse ot.txt (main output file with binding weights)
        ot_path = f"{self.working_dir}/ot.txt"
        if os.path.exists(ot_path):
            try:
                async with aiofiles.open(ot_path, 'r') as f:
                    ot_content = await f.read()
                    results["main_output"] = ot_content
                    
                    # Extract specific values
                    lines = ot_content.split('\n')
                    for i, line in enumerate(lines):
                        if "maxBindingWts" in line and i + 1 < len(lines):
                            weights_line = lines[i + 1].strip()
                            if weights_line:
                                try:
                                    weights = [float(x) for x in weights_line.split('\t') if x.strip()]
                                    results["binding_weights"] = weights
                                except ValueError:
                                    pass
                        elif "objective function" in line and i + 1 < len(lines):
                            obj_line = lines[i + 1].strip()
                            if obj_line:
                                try:
                                    results["objective_function"] = float(obj_line)
                                except ValueError:
                                    pass
            except Exception as e:
                results["ot_error"] = str(e)
        
        # Parse other output files if they exist
        for filename in ["ot3.txt", "pars2.txt"]:
            file_path = f"{self.working_dir}/{filename}"
            if os.path.exists(file_path):
                try:
                    async with aiofiles.open(file_path, 'r') as f:
                        content = await f.read()
                        if content.strip():  # Only include if not empty
                            results[filename] = content
                except Exception as e:
                    results[f"{filename}_error"] = str(e)
        
        return results
    
    async def _generate_pdf(self) -> Optional[str]:
        """Run scriptse.sh to convert format.tex to plot.pdf"""
        format_tex_path = f"{self.working_dir}/format.tex"
        script_path = SCRIPT_SE_PATH
        
        if not os.path.exists(format_tex_path):
            return None
            
        if not os.path.exists(script_path):
            return None
        
        try:
            # Run the script to generate PDF
            process = await asyncio.create_subprocess_exec(
                'bash', script_path,
                cwd=self.working_dir,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await asyncio.wait_for(
                process.communicate(), timeout=60  # 1 minute timeout for PDF generation
            )
            
            plot_pdf_path = f"{self.working_dir}/plot.pdf"
            if process.returncode == 0 and os.path.exists(plot_pdf_path):
                return plot_pdf_path
            else:
                return None
                
        except (asyncio.TimeoutError, Exception):
            return None


# FastAPI app
app = FastAPI(
    title="seq2exp MCP Server",
    description="Real seq2exp bioinformatics application wrapper with direct config file interface",
    version="1.0.0"
)


@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "name": "seq2exp MCP Server",
        "description": "Real bioinformatics application wrapper for gene expression prediction",
        "version": "1.0.0",
        "endpoints": {
            "predict": "/tools/seq2exp_predict",
            "validate": "/tools/seq2exp_validate_config", 
            "template": "/tools/seq2exp_config_template",
            "results": "/tools/seq2exp_get_results",
            "pdf": "/results/pdf/{filename}"
        }
    }


@app.post("/tools/seq2exp_predict")
async def seq2exp_predict(request: Seq2ExpRequest):
    """Main seq2exp prediction tool"""
    try:
        executor = Seq2ExpExecutor()
        result = await executor.run_prediction(request.config_content)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/tools/seq2exp_validate_config")
async def validate_config(request: Seq2ExpRequest):
    """Validate seq2exp configuration"""
    try:
        validation_result = Seq2ExpConfigManager.validate_config_content(request.config_content)
        return validation_result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/tools/seq2exp_config_template")
async def get_config_template():
    """Get default seq2exp configuration template"""
    try:
        template = await Seq2ExpConfigManager.get_default_config()
        
        return {
            "template": template,
            "description": "Default seq2exp configuration template. Modify parameters as needed.",
            "usage": "Copy this template and modify the parameters for your specific analysis"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/tools/seq2exp_get_results")
async def get_latest_results():
    """Get latest seq2exp results from output files"""
    try:
        executor = Seq2ExpExecutor()
        results = await executor._parse_seq2exp_output()
        
        # Check if plot.pdf exists
        plot_pdf_path = f"{WORKING_DIR}/plot.pdf"
        pdf_available = os.path.exists(plot_pdf_path)
        
        return {
            "results": results,
            "pdf_available": pdf_available,
            "pdf_url": "/results/pdf/plot.pdf" if pdf_available else None,
            "timestamp": time.time()
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/results/pdf/{filename}")
async def get_pdf(filename: str):
    """Serve generated PDF files"""
    pdf_path = f"{WORKING_DIR}/{filename}"
    if os.path.exists(pdf_path) and filename.endswith('.pdf'):
        return FileResponse(pdf_path, media_type="application/pdf", filename=filename)
    else:
        raise HTTPException(status_code=404, detail="PDF not found")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8083)