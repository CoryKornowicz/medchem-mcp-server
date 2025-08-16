# Medchem MCP Server Setup Guide

This guide will help you set up and integrate the Medchem MCP Server with Claude Desktop.

## What's Included

Your MCP server template includes:

- **3 Example Tools:**
  - `calculate_molecular_weight`: Calculate molecular weight from formula (e.g., "C6H12O6")
  - `get_molecule_info`: Get information about common molecules (aspirin, caffeine, glucose)
  - `echo`: Simple echo tool for testing

- **2 Resources:**
  - `medchem://molecules/aspirin`: Aspirin molecule data
  - `medchem://molecules/caffeine`: Caffeine molecule data

## Installation & Testing

The server has already been installed and tested successfully! You can run the test again anytime:

```bash
cd /Users/corykornowicz/Documents/LiteFold/medchem-mcp-server
python test_server.py
```

## Adding to Claude Desktop

### Step 1: Locate Claude Desktop Config

Find your Claude Desktop configuration file:

**On macOS:**
```bash
~/Library/Application Support/Claude/claude_desktop_config.json
```

**On Windows:**
```bash
%APPDATA%/Claude/claude_desktop_config.json
```

**On Linux:**
```bash
~/.config/Claude/claude_desktop_config.json
```

### Step 2: Add Your Server

Open the config file and add your server configuration. If the file doesn't exist, create it:

```json
{
  "mcpServers": {
    "medchem-mcp-server": {
      "command": "python",
      "args": ["-m", "medchem_mcp_server.server"],
      "cwd": "/Users/corykornowicz/Documents/LiteFold/medchem-mcp-server"
    }
  }
}
```

**Important:** Update the `cwd` path to match your actual project location if different.

### Step 3: Restart Claude Desktop

1. Completely quit Claude Desktop
2. Restart the application
3. The server should now be available

### Step 4: Test in Claude Desktop

Try these example prompts in Claude Desktop:

1. **Test the echo tool:**
   ```
   Use the echo tool to say "Hello from MCP!"
   ```

2. **Calculate molecular weight:**
   ```
   Calculate the molecular weight of glucose (C6H12O6)
   ```

3. **Get molecule information:**
   ```
   Tell me about aspirin using the molecule info tool
   ```

4. **List available resources:**
   ```
   What resources are available from the medchem server?
   ```

## Troubleshooting

### Server Not Appearing
- Ensure the `cwd` path in the config is correct
- Verify Python and dependencies are installed
- Check Claude Desktop logs for errors

### Test the Server Manually
```bash
cd /Users/corykornowicz/Documents/LiteFold/medchem-mcp-server
echo '{"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2024-11-05", "capabilities": {"roots": {"listChanged": true}}, "clientInfo": {"name": "test", "version": "1.0"}}}' | python -m medchem_mcp_server.server
```

### Check Dependencies
```bash
pip install -e .
```

## Extending the Server

### Adding New Tools

1. Add a new tool to the `handle_list_tools()` function
2. Implement the tool logic in `handle_call_tool()`
3. Test with `python test_server.py`

### Adding New Resources

1. Add resource metadata to `handle_list_resources()`
2. Implement resource reading in `handle_read_resource()`

### Example: Adding a New Tool

```python
# In handle_list_tools(), add:
Tool(
    name="convert_smiles_to_formula",
    description="Convert SMILES notation to molecular formula",
    inputSchema={
        "type": "object",
        "properties": {
            "smiles": {
                "type": "string",
                "description": "SMILES notation (e.g., 'CCO')",
            }
        },
        "required": ["smiles"],
    },
)

# In handle_call_tool(), add:
elif name == "convert_smiles_to_formula":
    smiles = arguments.get("smiles", "")
    # Your implementation here
    result = f"SMILES: {smiles}\nFormula: [calculated formula]"
    return [TextContent(type="text", text=result)]
```

## Next Steps

1. **Customize for your needs**: Replace the example tools with your specific medicinal chemistry tools
2. **Add real data**: Connect to chemical databases or calculation libraries
3. **Enhance functionality**: Add more sophisticated molecular calculations
4. **Error handling**: Improve error handling and validation

Your MCP server is now ready to use with Claude Desktop! 🎉
