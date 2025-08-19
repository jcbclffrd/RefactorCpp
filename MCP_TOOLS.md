# ExprPredictor MCP Tools

This directory contains MCP (Model Context Protocol) tools for the ExprPredictor library, enabling remote invocation of ExprPredictor functions through the GitHub MCP server.

## Overview

The MCP tools expose key functions from the ExprPredictor.cpp library as remotely callable tools. This allows agent workflows to interact with the gene expression prediction functionality without needing to compile or link against the C++ codebase directly.

## Available Tools

### 1. `expr_predictor_obj_func`
Computes the objective function value for given ExprPredictor parameters.

**Parameters:**
- `parameters`: ExprPar object containing model parameters

**Returns:**
- `success`: Boolean indicating success/failure
- `value`: Computed objective function value
- `message`: Status message

### 2. `expr_par_get_free_pars`
Extracts free parameters from an ExprPar object.

**Parameters:**
- `parameters`: ExprPar object
- `coopMat`: Cooperativity matrix (IntMatrix)
- `actIndicators`: Array of boolean activator indicators
- `repIndicators`: Array of boolean repressor indicators

**Returns:**
- `success`: Boolean indicating success/failure
- `freePars`: Array of extracted free parameters
- `message`: Status message

### 3. `expr_predictor_train`
Trains the ExprPredictor model with given data.

**Parameters:**
- `initialParameters`: Initial ExprPar parameters
- `trainingData`: Training dataset
- `options`: Training options

**Returns:**
- `success`: Boolean indicating success/failure
- `trainedParameters`: Trained model parameters
- `trainingStats`: Training statistics

### 4. `expr_predictor_predict`
Predicts expression values using a trained ExprPredictor.

**Parameters:**
- `parameters`: Trained ExprPar parameters
- `sequences`: Input sequences for prediction
- `conditions`: Experimental conditions

**Returns:**
- `success`: Boolean indicating success/failure
- `predictions`: Predicted expression values
- `message`: Status message

### 5. `expr_par_load`
Loads ExprPar parameters from a file.

**Parameters:**
- `filename`: Path to parameter file
- `coopMat`: Cooperativity matrix
- `repIndicators`: Repressor indicators

**Returns:**
- `success`: Boolean indicating success/failure
- `parameters`: Loaded ExprPar object
- `message`: Status message

### 6. `expr_func_predict_expr`
Predicts expression using ExprFunc.

**Parameters:**
- `sites`: Binding sites (SiteVec)
- `length`: Sequence length
- `factorConcentrations`: TF concentration values

**Returns:**
- `success`: Boolean indicating success/failure
- `expression`: Predicted expression value
- `message`: Status message

## Data Types

### ExprPar Object
```json
{
  "maxBindingWts": [array of doubles],
  "factorIntMat": {
    "rows": integer,
    "cols": integer,
    "data": [[array of arrays of doubles]]
  },
  "txpEffects": [array of doubles],
  "repEffects": [array of doubles],
  "basalTxp": double,
  "nFactors": integer
}
```

### Matrix Object
```json
{
  "rows": integer,
  "cols": integer,
  "data": [[array of arrays of doubles]]
}
```

### IntMatrix Object
```json
{
  "rows": integer,
  "cols": integer,
  "data": [[array of arrays of integers]]
}
```

## Building

```bash
make mcp_demo
```

This compiles the MCP tools along with the required ExprPredictor and Tools libraries.

## Testing

Run the demo program to test the MCP tools:

```bash
./mcp_demo
```

This will:
1. Initialize all MCP tools
2. List available tools
3. Show tool schemas
4. Execute sample tool calls

## Integration with GitHub MCP Server

The tools are designed to work with the GitHub MCP server for remote invocation. The MCP server can discover and call these tools through the standard MCP protocol.

### Example Tool Call

```json
{
  "method": "tools/call",
  "params": {
    "name": "expr_par_get_free_pars",
    "arguments": {
      "parameters": {
        "maxBindingWts": [1.0, 1.5, 2.0],
        "factorIntMat": {
          "rows": 3,
          "cols": 3,
          "data": [[1.0, 0.5, 0.0], [0.5, 1.0, 0.3], [0.0, 0.3, 1.0]]
        },
        "txpEffects": [1.2, 1.5, 0.8],
        "repEffects": [0.0, 0.0, 0.2],
        "basalTxp": 1.0
      },
      "coopMat": {
        "rows": 3,
        "cols": 3,
        "data": [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
      },
      "actIndicators": [true, true, false],
      "repIndicators": [false, false, true]
    }
  }
}
```

## Implementation Details

- **Language**: C++ with jsoncpp for JSON handling
- **Dependencies**: GSL (GNU Scientific Library), jsoncpp
- **Architecture**: Tool registry pattern with JSON serialization
- **Protocol**: MCP (Model Context Protocol) compatible

## Error Handling

All tools return a standardized response format with:
- `success`: Boolean indicating operation success
- `error`: Error message (if success is false)
- Additional fields specific to each tool

## Future Enhancements

1. Complete implementation of placeholder tools
2. Add more ExprPredictor functions as tools
3. Implement file I/O tools for loading/saving data
4. Add batch processing capabilities
5. Integrate with web server for HTTP access
6. Add parameter validation and type checking