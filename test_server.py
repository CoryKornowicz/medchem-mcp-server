#!/usr/bin/env python3
"""
Test script for the medchem MCP server.

This script tests the basic functionality of the MCP server to ensure it's working correctly.
"""

import asyncio
import json
import subprocess
import sys
from typing import Dict, Any


async def test_mcp_server():
    """Test the MCP server functionality."""
    print("🧪 Testing Medchem MCP Server...")
    print("=" * 50)
    
    # Test 1: Initialize and list tools
    print("\n1. Testing server initialization and tool listing...")
    
    init_request = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {
                "roots": {
                    "listChanged": True
                },
                "sampling": {}
            },
            "clientInfo": {
                "name": "test-client",
                "version": "1.0.0"
            }
        }
    }
    
    list_tools_request = {
        "jsonrpc": "2.0", 
        "id": 2,
        "method": "tools/list",
        "params": {}
    }
    
    # Notification to complete initialization
    initialized_notification = {
        "jsonrpc": "2.0",
        "method": "notifications/initialized",
        "params": {}
    }
    
    # Test 2: Call echo tool
    echo_request = {
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {
            "name": "echo",
            "arguments": {
                "message": "Hello from test!"
            }
        }
    }
    
    # Test 3: Calculate molecular weight
    molweight_request = {
        "jsonrpc": "2.0",
        "id": 4,
        "method": "tools/call",
        "params": {
            "name": "calculate_molecular_weight",
            "arguments": {
                "formula": "C6H12O6"
            }
        }
    }
    
    # Test 4: Get molecule info
    molinfo_request = {
        "jsonrpc": "2.0",
        "id": 5,
        "method": "tools/call",
        "params": {
            "name": "get_molecule_info",
            "arguments": {
                "name": "aspirin"
            }
        }
    }
    
    # Test 5: List resources
    list_resources_request = {
        "jsonrpc": "2.0",
        "id": 6,
        "method": "resources/list",
        "params": {}
    }
    
    try:
        # Start the server process
        print("Starting MCP server process...")
        process = await asyncio.create_subprocess_exec(
            sys.executable, "-m", "medchem_mcp_server.server",
            stdin=asyncio.subprocess.PIPE,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd="/Users/corykornowicz/Documents/LiteFold/medchem-mcp-server"
        )
        
        # Send requests and get responses
        requests = [
            init_request,
            initialized_notification,
            list_tools_request, 
            echo_request,
            molweight_request,
            molinfo_request,
            list_resources_request
        ]
        
        input_data = "\n".join(json.dumps(req) for req in requests) + "\n"
        
        print("Sending test requests...")
        stdout, stderr = await asyncio.wait_for(
            process.communicate(input_data.encode()),
            timeout=10.0
        )
        
        # Parse responses
        responses = []
        for line in stdout.decode().strip().split('\n'):
            if line.strip():
                try:
                    responses.append(json.loads(line))
                except json.JSONDecodeError as e:
                    print(f"Failed to parse response: {line}")
                    print(f"Error: {e}")
        
        # Analyze results
        print(f"\n✅ Received {len(responses)} responses")
        
        for i, response in enumerate(responses, 1):
            print(f"\nResponse {i}:")
            if "error" in response:
                print(f"  ❌ Error: {response['error']}")
            elif "result" in response:
                result = response["result"]
                if i == 1:  # Initialize response
                    print(f"  ✅ Server initialized: {result.get('serverInfo', {}).get('name', 'Unknown')}")
                elif i == 3:  # List tools (after init notification)
                    tools = result.get("tools", [])
                    print(f"  ✅ Found {len(tools)} tools:")
                    for tool in tools:
                        print(f"    - {tool['name']}: {tool['description']}")
                elif i == 4:  # Echo test
                    content = result.get("content", [])
                    if content and content[0].get("text"):
                        print(f"  ✅ Echo response: {content[0]['text']}")
                elif i == 5:  # Molecular weight
                    content = result.get("content", [])
                    if content and content[0].get("text"):
                        print(f"  ✅ Molecular weight calculation:")
                        print(f"    {content[0]['text']}")
                elif i == 6:  # Molecule info
                    content = result.get("content", [])
                    if content and content[0].get("text"):
                        print(f"  ✅ Molecule info:")
                        print(f"    {content[0]['text']}")
                elif i == 7:  # List resources
                    resources = result.get("resources", [])
                    print(f"  ✅ Found {len(resources)} resources:")
                    for resource in resources:
                        print(f"    - {resource['name']}")
            elif i == 2:  # Initialized notification (no response expected)
                print(f"  ✅ Initialization notification sent")
            else:
                print(f"  ⚠️  Unexpected response format: {response}")
        
        if stderr:
            print(f"\n⚠️  Server stderr: {stderr.decode()}")
        
        print(f"\n🎉 Test completed! Server is {'working correctly' if len(responses) >= 5 else 'partially working'}")
        
        return len(responses) >= 5
        
    except asyncio.TimeoutError:
        print("\n❌ Test timed out - server may be hanging")
        return False
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        return False
    finally:
        if 'process' in locals():
            try:
                process.terminate()
                await process.wait()
            except:
                pass


def main():
    """Run the test."""
    try:
        success = asyncio.run(test_mcp_server())
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\n❌ Test interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
