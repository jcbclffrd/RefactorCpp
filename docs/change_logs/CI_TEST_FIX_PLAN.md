# CI/CD Test Fix Plan

## Issue
GitHub Actions CI failing on `test_mcp_stdio.py` line 40 due to path resolution issues in the sandbox environment.

## Root Cause
Three tests use relative paths (`../seq2exp_mcp_server.py`) and `cwd=".."` which fails in GitHub Actions sandbox:
- `test_mcp_stdio.py` 
- `comprehensive_test_mcp.py`
- `test_validation_tool.py` (hardcoded path)

## Required Tests (from workflow)
1. `test_server_tools.py` ✅ (works - imports directly)
2. `test_mcp_stdio.py` ❌ (fails on line 40) 
3. `comprehensive_test_mcp.py` ❌ (same path issue)
4. `test_validation_tool.py` ❌ (hardcoded path issue)

## Solution
Replace relative paths with absolute path resolution:
```python
# Instead of: sys.executable, "../seq2exp_mcp_server.py", cwd=".."
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
server_script = os.path.join(parent_dir, "seq2exp_mcp_server.py")
# Use: sys.executable, server_script, cwd=parent_dir
```

## Archive Status
All other test files in `/Tests/` are legacy/archive - only the 4 above are actively used in CI.