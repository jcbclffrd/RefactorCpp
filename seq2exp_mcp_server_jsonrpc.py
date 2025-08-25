#!/usr/bin/env python3
"""
seq2exp MCP Server - True MCP protocol implementation
Uses JSON-RPC over stdio for Claude Desktop compatibility
"""

import asyncio
import json
import sys
import os
import tempfile
import time
from typing import Dict, Any, Optional, List

# Import your existing classes
from seq2exp_mcp_server import Seq2ExpExecutor, Seq2ExpConfigManager, WORKING_DIR


class MCPServer:
    """True MCP server implementation using JSON-RPC over stdio"""
    
    def __init__(self):
        self.executor = Seq2ExpExecutor()
        self.config_manager = Seq2ExpConfigManager()
    
    async def handle_message(self, message: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Handle incoming MCP JSON-RPC messages"""
        method = message.get("method")
        params = message.get("params", {})
        msg_id = message.get("id")
        
        try:
            if method == "initialize":
                return {
                    "jsonrpc": "2.0",
                    "id": msg_id,
                    "result": {
                        "protocolVersion": "2024-11-05",
                        "capabilities": {
                            "tools": {}
                        },
                        "serverInfo": {
                            "name": "seq2exp-mcp-server",
                            "version": "1.0.0"
                        }
                    }
                }
            
            elif method == "tools/list":
                return {
                    "jsonrpc": "2.0",
                    "id": msg_id,
                    "result": {
                        "tools": [
                            {
                                "name": "seq2exp_predict",
                                "description": "Run seq2exp prediction with configuration content",
                                "inputSchema": {
                                    "type": "object",
                                    "properties": {
                                        "config_content": {
                                            "type": "string",
                                            "description": "Complete seq2exp configuration file content"
                                        }
                                    },
                                    "required": ["config_content"]
                                }
                            },
                            {
                                "name": "seq2exp_validate_config",
                                "description": "Validate seq2exp configuration content",
                                "inputSchema": {
                                    "type": "object",
                                    "properties": {
                                        "config_content": {
                                            "type": "string",
                                            "description": "Configuration content to validate"
                                        }
                                    },
                                    "required": ["config_content"]
                                }
                            },
                            {
                                "name": "seq2exp_config_template",
                                "description": "Get default seq2exp configuration template",
                                "inputSchema": {
                                    "type": "object",
                                    "properties": {}
                                }
                            },
                            {
                                "name": "seq2exp_get_results",
                                "description": "Get latest seq2exp results from output files",
                                "inputSchema": {
                                    "type": "object",
                                    "properties": {}
                                }
                            }
                        ]
                    }
                }
            
            elif method == "tools/call":
                tool_name = params.get("name")
                arguments = params.get("arguments", {})
                
                if tool_name == "seq2exp_predict":
                    config_content = arguments.get("config_content", "")
                    result = await self.executor.run_prediction(config_content)
                    return {
                        "jsonrpc": "2.0",
                        "id": msg_id,
                        "result": {
                            "content": [
                                {
                                    "type": "text",
                                    "text": json.dumps(result, indent=2)
                                }
                            ]
                        }
                    }
                
                elif tool_name == "seq2exp_validate_config":
                    config_content = arguments.get("config_content", "")
                    validation_result = self.config_manager.validate_config_content(config_content)
                    return {
                        "jsonrpc": "2.0",
                        "id": msg_id,
                        "result": {
                            "content": [
                                {
                                    "type": "text",
                                    "text": json.dumps(validation_result, indent=2)
                                }
                            ]
                        }
                    }
                
                elif tool_name == "seq2exp_config_template":
                    template = await self.config_manager.get_default_config()
                    return {
                        "jsonrpc": "2.0",
                        "id": msg_id,
                        "result": {
                            "content": [
                                {
                                    "type": "text",
                                    "text": f"seq2exp Configuration Template:\n\n{template}"
                                }
                            ]
                        }
                    }
                
                elif tool_name == "seq2exp_get_results":
                    results = await self.executor._parse_seq2exp_output()
                    plot_pdf_path = f"{WORKING_DIR}/plot.pdf"
                    pdf_available = os.path.exists(plot_pdf_path)
                    
                    result_data = {
                        "results": results,
                        "pdf_available": pdf_available,
                        "timestamp": time.time()
                    }
                    
                    return {
                        "jsonrpc": "2.0",
                        "id": msg_id,
                        "result": {
                            "content": [
                                {
                                    "type": "text",
                                    "text": json.dumps(result_data, indent=2)
                                }
                            ]
                        }
                    }
                
                else:
                    raise ValueError(f"Unknown tool: {tool_name}")
            
            else:
                raise ValueError(f"Unknown method: {method}")
                
        except Exception as e:
            return {
                "jsonrpc": "2.0",
                "id": msg_id,
                "error": {
                    "code": -32603,
                    "message": str(e)
                }
            }
    
    async def run(self):
        """Main server loop - read from stdin, write to stdout"""
        while True:
            try:
                # Read line from stdin
                line = await asyncio.get_event_loop().run_in_executor(
                    None, sys.stdin.readline
                )
                
                if not line:
                    break
                
                # Parse JSON message
                try:
                    message = json.loads(line.strip())
                except json.JSONDecodeError:
                    continue
                
                # Handle message
                response = await self.handle_message(message)
                
                if response:
                    # Write response to stdout
                    print(json.dumps(response))
                    sys.stdout.flush()
                    
            except Exception as e:
                # Log error to stderr
                print(f"Error: {e}", file=sys.stderr)


if __name__ == "__main__":
    server = MCPServer()
    asyncio.run(server.run())
