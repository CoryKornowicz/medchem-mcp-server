
"""
Protein-related tools for the medchem MCP server.

This module contains all the protein structure and analysis tools.
"""

import logging
import requests
from typing import Optional
import mcp.types as types
from medchem_mcp_server.session_state import session_manager

logger = logging.getLogger("medchem-mcp-server")


def fetch_pdb_file(pdb_id: str, session_id: Optional[str] = None, 
                  save_to_session: bool = True, tags: Optional[str] = None) -> types.TextContent:
    """Fetch a PDB file from the RCSB Protein Data Bank and optionally save to session
    
    Args:
        pdb_id: The 4-character PDB identifier (e.g., '1abc')
        session_id: Optional session ID to save the file to (creates new if not provided)
        save_to_session: Whether to save the file to the session (default: True)
        tags: Optional comma-separated list of tags for the stored protein
    
    Returns:
        TextContent containing the PDB file contents or error message
    """
    # Clean and validate PDB ID
    pdb_id = pdb_id.strip().lower()

    if len(pdb_id) != 4:
        return types.TextContent(
            type="text",
            text="Error: PDB ID must be exactly 4 characters long"
        )

    # RCSB PDB download URL
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        return _fetch_pdb_file(pdb_url, pdb_id, session_id, save_to_session, tags)
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return types.TextContent(
                type="text",
                text=f"Error: PDB ID {pdb_id.upper()} not found in RCSB database"
            )
        else:
            return types.TextContent(
                type="text",
                text=f"Error: HTTP {e.response.status_code} when fetching PDB {pdb_id.upper()}"
            )
    except requests.exceptions.RequestException as e:
        return types.TextContent(
            type="text",
            text=f"Error: Network error when fetching PDB {pdb_id.upper()}: {str(e)}"
        )
    except Exception as e:
        logger.error(f"Unexpected error fetching PDB {pdb_id}: {e}", exc_info=True)
        return types.TextContent(
            type="text",
            text=f"Error: Unexpected error when fetching PDB {pdb_id.upper()}: {str(e)}"
        )


def _fetch_pdb_file(pdb_url, pdb_id, session_id=None, save_to_session=True, tags=None):
    response = requests.get(pdb_url, timeout=30)
    response.raise_for_status()

    # Check if we got a valid PDB file (should start with header info)
    pdb_content = response.text
    if not pdb_content.strip():
        return types.TextContent(
            type="text",
            text=f"Error: Empty response for PDB ID {pdb_id.upper()}"
        )

    # Basic validation - PDB files typically start with HEADER, TITLE, or similar
    first_line = pdb_content.split('\n')[0].strip()
    if not any(first_line.startswith(keyword) for keyword in ['HEADER', 'TITLE', 'REMARK', 'ATOM', 'HETATM']):
        return types.TextContent(
            type="text",
            text=f"Error: Invalid PDB file format for ID {pdb_id.upper()}"
        )

    # Save to session if requested
    session_info = ""
    if save_to_session:
        # Parse tags if provided
        tag_list = []
        if tags:
            tag_list = [tag.strip() for tag in tags.split(",")]
        
        # Get or create session
        session_id = session_manager.get_or_create_session(session_id)
        
        # Check if protein already exists in session
        existing = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
        if existing:
            protein_id, protein = existing
            session_info = f"\n\n[Protein already in session - Session ID: {session_id}, Protein ID: {protein_id[:8]}...]"
        else:
            # Save PDB file to temp directory
            file_id = session_manager.save_file(
                session_id=session_id,
                content=pdb_content,
                filename=f"{pdb_id}.pdb",
                file_type="pdb",
                tags=tag_list
            )
            
            file_entry = session_manager.get_file(session_id, file_id)
            file_path = file_entry.file_path if file_entry else None
            
            # Store protein in session
            protein_id = session_manager.add_protein(
                session_id=session_id,
                pdb_id=pdb_id,
                name=None,  # Could extract from TITLE record if needed
                pdb_content=pdb_content,
                file_path=file_path,
                tags=tag_list
            )
            
            # Get stored protein for details
            protein = session_manager.get_protein(session_id, protein_id)
            
            session_info = (
                f"\n\n[Saved to session - Session ID: {session_id}, "
                f"Protein ID: {protein_id[:8]}..., File ID: {file_id[:8]}..., "
                f"Chains: {protein.chain_count or 'N/A'}, "
                f"Residues: {protein.residue_count or 'N/A'}]"
            )

    return types.TextContent(
        type="text",
        text=f"PDB file for {pdb_id.upper()}:\n\n{session_info}"
    )
