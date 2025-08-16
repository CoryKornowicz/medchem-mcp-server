#!/usr/bin/env python3
"""
A simple MCP server template for medicinal chemistry applications.

This server provides basic tools for molecular calculations and data retrieval.
"""

import logging
import sys

from mcp.server.fastmcp import FastMCP

# Import all molecule tools
from medchem_mcp_server.molecule_tools import (
    get_molecule_descriptors_from_smiles,
    search_substructure,
    get_smiles_from_name,
    get_weight_from_smiles,
    find_maximum_common_substructure,
    smiles_to_image
)

# Import protein tools
from medchem_mcp_server.protein_tools import (
    fetch_pdb_file
)

# Import session-aware tools
from medchem_mcp_server.session_tools import (
    # Session management
    create_session,
    get_session_summary,
    clear_session,
    # Molecule management
    store_molecule,
    list_stored_molecules,
    get_stored_molecule,
    # Ligand conversion
    convert_ligand_pdb_to_sdf,
    convert_protein_ligands_to_sdf,
    get_ligand_image,
    # Docking preparation
    prepare_protein_for_docking,
    prepare_ligand_for_docking,
    # Protein management
    store_protein,
    list_stored_proteins,
    get_stored_protein,
    get_protein_component_files,
    analyze_stored_protein,
    # File management
    list_stored_files,
    get_stored_file,
    read_stored_file,
    # Collection management
    create_collection,
    get_collection
)

# Configure logging to stderr to avoid interfering with JSON-RPC over stdout
logging.basicConfig(
    level=logging.INFO,
    stream=sys.stderr,  # Send logs to stderr, not stdout
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger("medchem-mcp-server")

# Create FastMCP server instance
mcp = FastMCP("medchem-mcp-server")

# Register molecule tools
mcp.tool()(get_molecule_descriptors_from_smiles)
mcp.tool()(search_substructure)
mcp.tool()(get_smiles_from_name)
mcp.tool()(get_weight_from_smiles)
mcp.tool()(find_maximum_common_substructure)
mcp.tool()(smiles_to_image)

# Register protein tools
mcp.tool()(fetch_pdb_file)

# Register session management tools
mcp.tool()(create_session)
mcp.tool()(get_session_summary)
mcp.tool()(clear_session)

# Register molecule session tools
mcp.tool()(store_molecule)
mcp.tool()(list_stored_molecules)
mcp.tool()(get_stored_molecule)

# Register ligand conversion tools
mcp.tool()(convert_ligand_pdb_to_sdf)
mcp.tool()(convert_protein_ligands_to_sdf)
mcp.tool()(get_ligand_image)

# Register docking preparation tools
mcp.tool()(prepare_protein_for_docking)
mcp.tool()(prepare_ligand_for_docking)

# Register protein session tools
mcp.tool()(store_protein)
mcp.tool()(list_stored_proteins)
mcp.tool()(get_stored_protein)
mcp.tool()(get_protein_component_files)
mcp.tool()(analyze_stored_protein)

# Register file session tools
mcp.tool()(list_stored_files)
mcp.tool()(get_stored_file)
mcp.tool()(read_stored_file)

# Register collection tools
mcp.tool()(create_collection)
mcp.tool()(get_collection)


if __name__ == "__main__":
    try:
        mcp.run()
    except Exception as e:
        logger.error(f"Server error: {e}", exc_info=True)
        sys.exit(1)