#!/usr/bin/env python3
"""
A simple MCP server template for medicinal chemistry applications.

This server provides basic tools for molecular calculations and data retrieval.
"""

import logging
import re
from typing import Any

from mcp.server.fastmcp import FastMCP

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("medchem-mcp-server")

# Create FastMCP server instance
mcp = FastMCP("medchem-mcp-server")


# Resource data
MOLECULE_DATA = {
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
    }
}

@mcp.resource("medchem://molecules/{molecule}")
def get_molecule_resource(molecule: str) -> str:
    """Get molecule data as a JSON resource."""
    if molecule not in MOLECULE_DATA:
        raise ValueError(f"Unknown molecule: {molecule}")
    
    import json
    return json.dumps(MOLECULE_DATA[molecule], indent=2)


# Atomic weights for molecular weight calculations
ATOMIC_WEIGHTS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904,
    'F': 18.998, 'I': 126.904, 'Na': 22.990, 'K': 39.098,
    'Ca': 40.078, 'Mg': 24.305, 'Fe': 55.845, 'Zn': 65.38,
}


@mcp.tool()
def calculate_molecular_weight(formula: str) -> str:
    """Calculate molecular weight from a molecular formula.
    
    Args:
        formula: Molecular formula (e.g., 'C6H12O6')
        
    Returns:
        Detailed breakdown of molecular weight calculation
    """
    if not formula:
        raise ValueError("Formula is required")
        
    # Basic parser for simple formulas like C6H12O6
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    
    total_weight = 0.0
    composition = []
    
    for element, count in matches:
        count = int(count) if count else 1
        if element in ATOMIC_WEIGHTS:
            weight = ATOMIC_WEIGHTS[element] * count
            total_weight += weight
            composition.append(f"{element}: {count} atoms, {weight:.3f} g/mol")
        else:
            raise ValueError(f"Unknown element: {element}")
    
    result = f"Molecular formula: {formula}\n"
    result += f"Molecular weight: {total_weight:.3f} g/mol\n\n"
    result += "Composition:\n" + "\n".join(composition)
    
    return result


@mcp.tool()
def get_molecule_info(name: str) -> str:
    """Get basic information about a molecule by name.
    
    Args:
        name: Name of the molecule (e.g., 'aspirin', 'caffeine')
        
    Returns:
        Detailed molecule information including formula, weight, SMILES, and description
    """
    molecule_name = name.lower()
    
    # Extended molecule database
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
        return result
    else:
        available = ", ".join(molecules.keys())
        return f"Molecule '{molecule_name}' not found. Available molecules: {available}"


@mcp.tool()
def echo(message: str) -> str:
    """Simple echo tool for testing.
    
    Args:
        message: Message to echo back
        
    Returns:
        The echoed message with prefix
    """
    return f"Echo: {message}"


if __name__ == "__main__":
    mcp.run()
