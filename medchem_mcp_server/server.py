#!/usr/bin/env python3
"""
A simple MCP server template for medicinal chemistry applications.

This server provides basic tools for molecular calculations and data retrieval.
"""

import asyncio
import logging
from typing import Any, Sequence

from mcp.server.models import InitializationOptions
from mcp.server import NotificationOptions, Server
from mcp.types import (
    Resource,
    Tool,
    TextContent,
    ImageContent,
    EmbeddedResource,
    LoggingLevel
)
from pydantic import AnyUrl
import mcp.server.stdio


# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("medchem-mcp-server")

# Create server instance
server = Server("medchem-mcp-server")


@server.list_resources()
async def handle_list_resources() -> list[Resource]:
    """List available resources."""
    return [
        Resource(
            uri=AnyUrl("medchem://molecules/aspirin"),
            name="Aspirin molecule data",
            description="Basic information about aspirin",
            mimeType="application/json",
        ),
        Resource(
            uri=AnyUrl("medchem://molecules/caffeine"),
            name="Caffeine molecule data", 
            description="Basic information about caffeine",
            mimeType="application/json",
        ),
    ]


@server.read_resource()
async def handle_read_resource(uri: AnyUrl) -> str:
    """Read a specific resource."""
    if uri.scheme != "medchem":
        raise ValueError(f"Unsupported URI scheme: {uri.scheme}")
    
    path = str(uri).replace("medchem://", "")
    
    if path == "molecules/aspirin":
        return """
        {
            "name": "Aspirin",
            "formula": "C9H8O4",
            "molecular_weight": 180.16,
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "description": "Aspirin is a medication used to reduce pain, fever, or inflammation."
        }
        """
    elif path == "molecules/caffeine":
        return """
        {
            "name": "Caffeine",
            "formula": "C8H10N4O2", 
            "molecular_weight": 194.19,
            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "description": "Caffeine is a central nervous system stimulant."
        }
        """
    else:
        raise ValueError(f"Unknown resource path: {path}")


@server.list_tools()
async def handle_list_tools() -> list[Tool]:
    """List available tools."""
    return [
        Tool(
            name="calculate_molecular_weight",
            description="Calculate molecular weight from a molecular formula",
            inputSchema={
                "type": "object",
                "properties": {
                    "formula": {
                        "type": "string",
                        "description": "Molecular formula (e.g., 'C6H12O6')",
                    }
                },
                "required": ["formula"],
            },
        ),
        Tool(
            name="get_molecule_info",
            description="Get basic information about a molecule by name",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Name of the molecule (e.g., 'aspirin', 'caffeine')",
                    }
                },
                "required": ["name"],
            },
        ),
        Tool(
            name="echo",
            description="Simple echo tool for testing",
            inputSchema={
                "type": "object",
                "properties": {
                    "message": {
                        "type": "string",
                        "description": "Message to echo back",
                    }
                },
                "required": ["message"],
            },
        ),
    ]


@server.call_tool()
async def handle_call_tool(name: str, arguments: dict[str, Any] | None) -> list[TextContent | ImageContent | EmbeddedResource]:
    """Handle tool calls."""
    if arguments is None:
        arguments = {}
        
    if name == "calculate_molecular_weight":
        formula = arguments.get("formula", "")
        if not formula:
            raise ValueError("Formula is required")
            
        # Simple molecular weight calculation (basic implementation)
        atomic_weights = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904,
            'F': 18.998, 'I': 126.904, 'Na': 22.990, 'K': 39.098,
            'Ca': 40.078, 'Mg': 24.305, 'Fe': 55.845, 'Zn': 65.38,
        }
        
        # Basic parser for simple formulas like C6H12O6
        import re
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        total_weight = 0.0
        composition = []
        
        for element, count in matches:
            count = int(count) if count else 1
            if element in atomic_weights:
                weight = atomic_weights[element] * count
                total_weight += weight
                composition.append(f"{element}: {count} atoms, {weight:.3f} g/mol")
            else:
                raise ValueError(f"Unknown element: {element}")
        
        result = f"Molecular formula: {formula}\n"
        result += f"Molecular weight: {total_weight:.3f} g/mol\n\n"
        result += "Composition:\n" + "\n".join(composition)
        
        return [TextContent(type="text", text=result)]
        
    elif name == "get_molecule_info":
        molecule_name = arguments.get("name", "").lower()
        
        molecules = {
            "aspirin": {
                "name": "Aspirin",
                "formula": "C9H8O4",
                "molecular_weight": 180.16,
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "description": "Aspirin is a medication used to reduce pain, fever, or inflammation."
            },
            "caffeine": {
                "name": "Caffeine",
                "formula": "C8H10N4O2",
                "molecular_weight": 194.19,
                "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "description": "Caffeine is a central nervous system stimulant."
            },
            "glucose": {
                "name": "Glucose",
                "formula": "C6H12O6",
                "molecular_weight": 180.16,
                "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
                "description": "Glucose is a simple sugar and an important carbohydrate in biology."
            }
        }
        
        if molecule_name in molecules:
            mol = molecules[molecule_name]
            result = f"Molecule: {mol['name']}\n"
            result += f"Formula: {mol['formula']}\n"
            result += f"Molecular Weight: {mol['molecular_weight']} g/mol\n"
            result += f"SMILES: {mol['smiles']}\n"
            result += f"Description: {mol['description']}"
            return [TextContent(type="text", text=result)]
        else:
            available = ", ".join(molecules.keys())
            return [TextContent(type="text", text=f"Molecule '{molecule_name}' not found. Available molecules: {available}")]
            
    elif name == "echo":
        message = arguments.get("message", "")
        return [TextContent(type="text", text=f"Echo: {message}")]
        
    else:
        raise ValueError(f"Unknown tool: {name}")


async def main():
    """Main entry point for the server."""
    # Run the server using stdin/stdout streams
    async with mcp.server.stdio.stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="medchem-mcp-server",
                server_version="0.1.0",
                capabilities=server.get_capabilities(
                    notification_options=NotificationOptions(),
                    experimental_capabilities={},
                ),
            ),
        )


if __name__ == "__main__":
    asyncio.run(main())
