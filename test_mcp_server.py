#!/usr/bin/env python3
"""
Test script for seq2exp MCP server
"""

import asyncio
import json
import time
import aiohttp

SERVER_URL = "http://localhost:8083"

async def test_server():
    """Test all MCP server endpoints"""
    
    print("Testing seq2exp MCP Server")
    print("=" * 50)
    
    async with aiohttp.ClientSession() as session:
        
        # Test 1: Get config template
        print("\n1. Testing config template endpoint...")
        try:
            async with session.get(f"{SERVER_URL}/tools/seq2exp_config_template") as resp:
                if resp.status == 200:
                    data = await resp.json()
                    print("✓ Config template retrieved successfully")
                    print(f"  Template length: {len(data['template'])} characters")
                else:
                    print(f"✗ Config template failed: {resp.status}")
                    print(await resp.text())
        except Exception as e:
            print(f"✗ Config template error: {e}")
        
        # Test 2: Validate config
        print("\n2. Testing config validation...")
        test_config = """# Test config
seqFile = iData/rhoseq.txt
exprFile = iData/rhoexp.tab
motifFile = iData/factordts.wtmx
factorExprFile = iData/factorexpdts2s.tab
outputFile = out.txt
nRandomStarts = 1
energyThreshold = 5"""
        
        try:
            async with session.post(f"{SERVER_URL}/tools/seq2exp_validate_config", 
                                  json={"config_content": test_config}) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    print(f"✓ Config validation completed")
                    print(f"  Valid: {data['valid']}")
                    print(f"  Issues: {len(data['issues'])}")
                    if data['issues']:
                        for issue in data['issues'][:3]:  # Show first 3 issues
                            print(f"    - {issue}")
                else:
                    print(f"✗ Config validation failed: {resp.status}")
                    print(await resp.text())
        except Exception as e:
            print(f"✗ Config validation error: {e}")
        
        # Test 3: Run prediction with minimal config
        print("\n3. Testing seq2exp prediction (minimal config)...")
        minimal_config = """# Fast test config for seq2exp
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
        
        try:
            start_time = time.time()
            async with session.post(f"{SERVER_URL}/tools/seq2exp_predict",
                                  json={"config_content": minimal_config}) as resp:
                request_time = time.time() - start_time
                
                if resp.status == 200:
                    data = await resp.json()
                    print(f"✓ Prediction completed in {request_time:.2f}s")
                    print(f"  Success: {data['success']}")
                    print(f"  Execution time: {data.get('execution_time', 'N/A'):.2f}s")
                    if data['success']:
                        results = data.get('results', {})
                        print(f"  Results available: {len(results)} items")
                        if 'binding_weights' in results:
                            print(f"  Binding weights: {results['binding_weights'][:3]}...")
                        if 'objective_function' in results:
                            print(f"  Objective function: {results['objective_function']}")
                        if data.get('pdf_path'):
                            print(f"  PDF generated: {data['pdf_path']}")
                    else:
                        print(f"  Error: {data.get('error', 'Unknown error')}")
                else:
                    print(f"✗ Prediction failed: {resp.status}")
                    error_text = await resp.text()
                    print(f"  Error: {error_text[:200]}...")
        except Exception as e:
            print(f"✗ Prediction error: {e}")
        
        # Test 4: Get latest results
        print("\n4. Testing get results endpoint...")
        try:
            async with session.get(f"{SERVER_URL}/tools/seq2exp_get_results") as resp:
                if resp.status == 200:
                    data = await resp.json()
                    print("✓ Results retrieved successfully")
                    print(f"  PDF available: {data['pdf_available']}")
                    if data['pdf_available']:
                        print(f"  PDF URL: {data['pdf_url']}")
                    results = data.get('results', {})
                    print(f"  Result items: {len(results)}")
                else:
                    print(f"✗ Get results failed: {resp.status}")
                    print(await resp.text())
        except Exception as e:
            print(f"✗ Get results error: {e}")
        
        # Test 5: Check if PDF endpoint works
        print("\n5. Testing PDF endpoint...")
        try:
            async with session.get(f"{SERVER_URL}/results/pdf/plot.pdf") as resp:
                if resp.status == 200:
                    content_length = resp.headers.get('content-length', 'unknown')
                    print(f"✓ PDF accessible (size: {content_length} bytes)")
                elif resp.status == 404:
                    print("- PDF not found (expected if no recent run)")
                else:
                    print(f"✗ PDF endpoint failed: {resp.status}")
        except Exception as e:
            print(f"✗ PDF endpoint error: {e}")

if __name__ == "__main__":
    asyncio.run(test_server())