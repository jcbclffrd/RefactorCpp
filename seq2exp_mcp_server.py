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
from fastmcp import FastMCP


# Configuration
# Dynamically determine the working directory
WORKING_DIR = os.getenv("WORKING_DIR", os.getcwd())
SEQ2EXP_EXECUTABLE = f"{WORKING_DIR}/seq2exp"
SCRIPT_SE_PATH = f"{WORKING_DIR}/scriptse.sh"
DEFAULT_CONFIG_PATH = f"{WORKING_DIR}/seq2exp.conf"


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
                # PDF generation failed, but don't treat as fatal error
                return None
                
        except (asyncio.TimeoutError, Exception):
            # PDF generation failed, but don't treat as fatal error
            return None


# FastMCP server
mcp = FastMCP(
    name="seq2exp MCP Server",
    version="1.0.0"
)


async def seq2exp_predict_impl(config_content: str) -> Dict[str, Any]:
    """
    Main seq2exp prediction tool.
    
    Args:
        config_content: Complete seq2exp configuration file content
        
    Returns:
        Dictionary containing prediction results, execution time, and output files
    """
    try:
        executor = Seq2ExpExecutor()
        result = await executor.run_prediction(config_content)
        return result
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }


async def seq2exp_validate_config_impl(config_content: str) -> Dict[str, Any]:
    """
    Validate seq2exp configuration.
    
    Args:
        config_content: seq2exp configuration file content to validate
        
    Returns:
        Dictionary containing validation results and any issues found
    """
    try:
        validation_result = Seq2ExpConfigManager.validate_config_content(config_content)
        return validation_result
    except Exception as e:
        return {
            "valid": False,
            "error": str(e)
        }


async def seq2exp_config_template_impl() -> Dict[str, Any]:
    """
    Get default seq2exp configuration template.
    
    Returns:
        Dictionary containing the default configuration template and usage instructions
    """
    try:
        template = await Seq2ExpConfigManager.get_default_config()
        
        return {
            "template": template,
            "description": "Default seq2exp configuration template. Modify parameters as needed.",
            "usage": "Copy this template and modify the parameters for your specific analysis"
        }
    except Exception as e:
        return {
            "error": str(e)
        }


async def seq2exp_get_results_impl() -> Dict[str, Any]:
    """
    Get latest seq2exp results from output files.
    
    Returns:
        Dictionary containing the latest results and PDF availability information
    """
    try:
        executor = Seq2ExpExecutor()
        results = await executor._parse_seq2exp_output()
        
        # Check if plot.pdf exists
        plot_pdf_path = f"{WORKING_DIR}/plot.pdf"
        pdf_available = os.path.exists(plot_pdf_path)
        
        return {
            "results": results,
            "pdf_available": pdf_available,
            "pdf_path": plot_pdf_path if pdf_available else None,
            "timestamp": time.time()
        }
    except Exception as e:
        return {
            "error": str(e)
        }


# Register the functions as MCP tools
mcp.tool(seq2exp_predict_impl, name="seq2exp_predict")
mcp.tool(seq2exp_validate_config_impl, name="seq2exp_validate_config") 
mcp.tool(seq2exp_config_template_impl, name="seq2exp_config_template")
mcp.tool(seq2exp_get_results_impl, name="seq2exp_get_results")


if __name__ == "__main__":
    # Run the MCP server using stdio transport for Claude Desktop compatibility
    asyncio.run(mcp.run_stdio_async())