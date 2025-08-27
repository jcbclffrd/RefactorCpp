# seq2exp MCP Server Setup Guide

This guide explains how to set up and use the seq2exp MCP Server with various hosts including Claude Desktop, VS Code, and Cursor.

## Overview

The seq2exp MCP Server provides bioinformatics tools through the Model Context Protocol (MCP). It wraps the seq2exp application and provides tools for:

- `seq2exp_predict_impl`: Run seq2exp predictions with configuration
- `seq2exp_validate_config_impl`: Validate seq2exp configuration files
- `seq2exp_config_template_impl`: Get configuration templates
- `seq2exp_get_results_impl`: Retrieve latest results and PDFs

## Prerequisites

- seq2exp MCP Server dependencies installed:
  ```bash
  uv pip install -r requirements_mcp.txt
  ```
- seq2exp executable and required data files in the project directory

## Claude Desktop Setup

### Windows with WSL Configuration

If you're running Claude Desktop on Windows but have the MCP server in WSL, use this configuration:

**File**: `C:\Users\[username]\AppData\Roaming\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "seq2exp": {
      "command": "wsl",
      "args": [
        "env",
        "-i",
        "bash",
        "-c",
        "cd /home/j/MPA && /home/j/miniconda/bin/python3 seq2exp_mcp_server.py"
      ]
    }
  }
}
```

**Key points:**
- `env -i` creates a clean environment, avoiding Windows PATH conflicts
- Full path to python3 ensures correct Python interpreter
- Adjust paths (`/home/j/MPA`, `/home/j/miniconda/bin/python3`) to match your setup

### Linux/macOS Configuration

**File**: `~/.config/claude-desktop/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "seq2exp": {
      "command": "python3",
      "args": ["/path/to/seq2exp_mcp_server.py"],
      "cwd": "/path/to/MPA"
    }
  }
}
```

### After Configuration

1. Save the configuration file
2. Completely close Claude Desktop
3. Restart Claude Desktop
4. The MCP server tools should now be available in new chats

## VS Code Setup

### Using MCP Extension

1. Install the MCP extension for VS Code
2. Configure in VS Code settings or `.vscode/settings.json`:

```json
{
  "mcp.servers": {
    "seq2exp": {
      "command": "python3",
      "args": ["/path/to/seq2exp_mcp_server.py"],
      "cwd": "/path/to/MPA"
    }
  }
}
```

### Using Claude Dev Extension

If using the Claude Dev extension:

1. Configure the MCP server in the extension settings
2. Or add to your workspace `.vscode/settings.json`:

```json
{
  "claude-dev.mcpServers": {
    "seq2exp": {
      "command": "python3",
      "args": ["/path/to/seq2exp_mcp_server.py"],
      "cwd": "/path/to/MPA"
    }
  }
}
```

## Cursor Setup

### Method 1: Direct MCP Configuration

Configure in Cursor's settings or `.cursor/settings.json`:

```json
{
  "mcp.servers": {
    "seq2exp": {
      "command": "python3",
      "args": ["/path/to/seq2exp_mcp_server.py"],
      "cwd": "/path/to/MPA"
    }
  }
}
```

### Method 2: Using Claude Integration

If Cursor supports Claude's MCP format, use the same configuration as Claude Desktop but in Cursor's config location.

## Other Host Applications

For any application that supports MCP servers, use this general pattern:

```json
{
  "command": "python3",
  "args": ["/path/to/seq2exp_mcp_server.py"],
  "cwd": "/path/to/MPA",
  "env": {
    "WORKING_DIR": "/path/to/MPA"
  }
}
```

## Troubleshooting

### Common Issues

1. **Module not found errors**:
   - Ensure dependencies are installed: `uv pip install -r requirements_mcp.txt`
   - Use full path to Python interpreter
   - For WSL: Use `env -i` to avoid PATH conflicts

2. **Server fails to start**:
   - Check that `seq2exp_mcp_server.py` is executable: `chmod +x seq2exp_mcp_server.py`
   - Verify Python path and working directory
   - Check application logs for specific errors

3. **WSL/Windows path issues**:
   - Use `env -i` to start with clean environment
   - Use full WSL paths, not Windows paths
   - Ensure WSL can access all required files

### Testing the Server

Test the server directly:

```bash
# Direct test
python3 seq2exp_mcp_server.py

# WSL test (from Windows)
wsl env -i bash -c "cd /home/j/MPA && /home/j/miniconda/bin/python3 seq2exp_mcp_server.py"
```

The server should start and show the FastMCP banner with "Starting MCP server" message.

### Logging

Most MCP clients provide logging. Check:
- Claude Desktop: Look for MCP server logs in the application logs
- VS Code: Check the Output panel for MCP-related messages
- Cursor: Check developer tools or output panels

## Usage Examples

Once connected, you can use the tools in natural language:

- "Validate this seq2exp configuration: [paste config]"
- "Run a seq2exp prediction with these parameters..."
- "Get the default seq2exp configuration template"
- "Show me the latest seq2exp results"

## Environment Variables

The server recognizes these environment variables:

- `WORKING_DIR`: Override the working directory (default: current directory)
- Additional variables can be set in the MCP configuration's `env` section

## Dependencies

The server requires these Python packages (see `requirements_mcp.txt`):
- `fastmcp>=2.11.0`
- `aiofiles>=23.2.0`
- `pydantic>=2.0.0,<3.0.0`

## File Structure

Expected project structure:
```
MPA/
├── seq2exp_mcp_server.py      # MCP server script
├── seq2exp                    # seq2exp executable
├── seq2exp.conf               # Default configuration
├── scriptse.sh                # PDF generation script
├── requirements_mcp.txt       # Python dependencies
└── iData/                     # Input data directory
```

## Security Notes

- The server executes seq2exp with user-provided configurations
- Ensure the working directory has appropriate permissions
- Consider running in a sandboxed environment for untrusted inputs
- File paths are validated but users can read/write within the working directory