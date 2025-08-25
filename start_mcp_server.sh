#!/bin/bash

# seq2exp MCP Server Startup Script
# This script sets up and starts the seq2exp MCP server

set -e

echo "seq2exp MCP Server Setup"
echo "========================"

# Check if we're in the right directory
if [ ! -f "seq2exp" ]; then
    echo "Error: seq2exp executable not found. Please run from the MPA repository root."
    exit 1
fi

if [ ! -d "iData" ]; then
    echo "Error: iData directory not found. Please run from the MPA repository root."
    exit 1
fi

# Check if seq2exp works
echo "Testing seq2exp executable..."
if ! ./seq2exp >/dev/null 2>&1; then
    echo "Warning: seq2exp may need GSL library. Installing..."
    sudo apt-get update && sudo apt-get install -y libgsl-dev
    
    echo "Rebuilding seq2exp..."
    make clean && make
fi

# Check Python dependencies
echo "Checking Python dependencies..."
python3 -c "import fastapi, uvicorn, aiofiles" 2>/dev/null || {
    echo "Installing Python dependencies..."
    python3 -m pip3 install fastapi uvicorn aiofiles python-multipart aiohttp
}

# Check LaTeX for PDF generation
echo "Checking LaTeX installation..."
if ! command -v latex >/dev/null 2>&1; then
    echo "Installing LaTeX for PDF generation..."
    sudo apt-get install -y texlive texlive-extra-utils
fi

echo ""
echo "Starting seq2exp MCP Server..."
echo "Server will be available at: http://localhost:8083"
echo "Interactive documentation: http://localhost:8083/docs"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Start the server
python3 seq2exp_mcp_server.py