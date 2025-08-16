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
    # Docking tools
    dock_ligand_with_vina,
    score_ligand_with_vina,
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
    get_collection,
    # Molecule editing
    edit_stored_molecule,
    template_stored_molecule,
    get_available_r_groups,
    get_available_scaffolds
)

# Import ZINC tools
from medchem_mcp_server.zinc_tooling import (
    fetch_zinc_molecule,
    fetch_zinc_batch,
    search_zinc_similarity,
    search_zinc_substructure,
    prepare_zinc_molecule_for_editing
)

# Import ZINC editing integration tools
from medchem_mcp_server.zinc_editing_tools import (
    fetch_and_edit_zinc_molecule,
    create_zinc_analog_series,
    scaffold_hop_zinc_molecule,
    fragment_and_recombine_zinc,
    apply_medchem_transforms,
    visualize_zinc_edits
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

# Register docking tools
mcp.tool()(dock_ligand_with_vina)
mcp.tool()(score_ligand_with_vina)

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

# Register molecule editing session tools
mcp.tool()(edit_stored_molecule)
mcp.tool()(template_stored_molecule)
mcp.tool()(get_available_r_groups)
mcp.tool()(get_available_scaffolds)

# Register ZINC fetching tools
mcp.tool()(fetch_zinc_molecule)
mcp.tool()(fetch_zinc_batch)
mcp.tool()(search_zinc_similarity)
mcp.tool()(search_zinc_substructure)
mcp.tool()(prepare_zinc_molecule_for_editing)

# Register ZINC editing integration tools
mcp.tool()(fetch_and_edit_zinc_molecule)
mcp.tool()(create_zinc_analog_series)
mcp.tool()(scaffold_hop_zinc_molecule)
mcp.tool()(fragment_and_recombine_zinc)
mcp.tool()(apply_medchem_transforms)
mcp.tool()(visualize_zinc_edits)


if __name__ == "__main__":
    try:
        mcp.run()
    except Exception as e:
        logger.error(f"Server error: {e}", exc_info=True)
        sys.exit(1)