"""
Session-aware tools for molecule and protein management.
"""

import logging
import os
from typing import Optional, List, Union
import mcp.types as types
from medchem_mcp_server.session_state import session_manager
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import json

# Import protein analysis utilities (lazy import to avoid dependency issues)
UTILITIES_AVAILABLE = False
try:
    from medchem_mcp_server.utilities import (
        get_protein,
        get_ligands,
        get_waters,
        get_ions,
        get_nucleic_acids,
        get_non_standard_amino_acids
    )
    UTILITIES_AVAILABLE = True
except ImportError:
    # Utilities not available, will use basic analysis
    pass

# Import molecule utilities
from medchem_mcp_server.molecule_utilities import (
    pdb_to_sdf,
    pdb_file_to_sdf_file,
    extract_smiles_from_pdb,
    extract_smiles_from_sdf,
    calculate_ligand_properties,
    ligand_atomgroup_to_sdf,
    process_ligand_stereochemistry,
    write_protein_pdbqt,
    write_ligand_pdbqt,
    add_missing_protein_atoms
)

logger = logging.getLogger("medchem-mcp-server")

# Molecule management tools

def store_molecule(smiles: str, name: Optional[str] = None, session_id: Optional[str] = None, 
                  tags: Optional[str] = None) -> types.TextContent:
    """Store a molecule in the session state for later use
    
    Args:
        smiles: SMILES string of the molecule
        name: Optional name for the molecule
        session_id: Optional session ID (will create new if not provided)
        tags: Comma-separated list of tags
    
    Returns:
        TextContent with storage confirmation and IDs
    """
    
    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return types.TextContent(
            type="text",
            text=f"Error: Invalid SMILES string: {smiles}"
        )
    
    # Parse tags if provided
    tag_list = []
    if tags:
        tag_list = [tag.strip() for tag in tags.split(",")]
    
    # Get or create session
    session_id = session_manager.get_or_create_session(session_id)
    
    # Calculate basic properties
    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": Descriptors.TPSA(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds()
    }
    
    # Store molecule
    molecule_id = session_manager.add_molecule(
        session_id=session_id,
        smiles=smiles,
        name=name,
        properties=properties,
        tags=tag_list
    )
    
    return types.TextContent(
        type="text",
        text=f"Molecule stored successfully!\nSession ID: {session_id}\nMolecule ID: {molecule_id}\nSMILES: {smiles}\nName: {name or 'N/A'}\nTags: {', '.join(tag_list) if tag_list else 'None'}"
    )

def list_stored_molecules(session_id: str, tags: Optional[str] = None) -> types.TextContent:
    """List molecules stored in the session
    
    Args:
        session_id: Session ID to query
        tags: Optional comma-separated list of tags to filter by
    
    Returns:
        TextContent with list of stored molecules
    """
    
    # Parse tags if provided
    tag_filter = None
    if tags:
        tag_filter = [tag.strip() for tag in tags.split(",")]
    
    molecules = session_manager.list_molecules(session_id, tag_filter)
    
    if not molecules:
        return types.TextContent(
            type="text",
            text="No molecules found in this session."
        )
    
    # Create summary table
    lines = ["Stored Molecules:", "=" * 50]
    for mol_id, molecule in molecules.items():
        lines.append(f"ID: {mol_id[:8]}...")
        lines.append(f"  Name: {molecule.name or 'N/A'}")
        lines.append(f"  SMILES: {molecule.smiles}")
        if 'molecular_weight' in molecule.properties:
            lines.append(f"  MW: {molecule.properties['molecular_weight']:.2f}")
        if 'logp' in molecule.properties:
            lines.append(f"  LogP: {molecule.properties['logp']:.2f}")
        lines.append(f"  Tags: {', '.join(molecule.tags) if molecule.tags else 'None'}")
        lines.append(f"  Created: {molecule.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("-" * 30)
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

def get_stored_molecule(session_id: str, molecule_id: str) -> types.TextContent:
    """Get details of a specific stored molecule
    
    Args:
        session_id: Session ID
        molecule_id: Molecule ID to retrieve
    
    Returns:
        TextContent with molecule details
    """
    
    molecule = session_manager.get_molecule(session_id, molecule_id)
    
    if not molecule:
        return types.TextContent(
            type="text",
            text=f"Molecule with ID {molecule_id} not found in session {session_id}"
        )
    
    details = [
        f"Molecule Details:",
        f"ID: {molecule_id}",
        f"Name: {molecule.name or 'N/A'}",
        f"SMILES: {molecule.smiles}",
        f"Tags: {', '.join(molecule.tags) if molecule.tags else 'None'}",
        f"Created: {molecule.created_at.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "Properties:",
    ]
    
    for prop, value in molecule.properties.items():
        if isinstance(value, float):
            details.append(f"  {prop}: {value:.3f}")
        else:
            details.append(f"  {prop}: {value}")
    
    return types.TextContent(
        type="text",
        text="\n".join(details)
    )

# Protein management tools

def store_protein(pdb_id: str, name: Optional[str] = None, session_id: Optional[str] = None,
                 pdb_content: Optional[str] = None, tags: Optional[str] = None) -> types.TextContent:
    """Store a protein in the session state
    
    Args:
        pdb_id: PDB identifier
        name: Optional name for the protein
        session_id: Optional session ID (will create new if not provided)
        pdb_content: Optional PDB file content
        tags: Comma-separated list of tags
    
    Returns:
        TextContent with storage confirmation and IDs
    """
    
    # Parse tags if provided
    tag_list = []
    if tags:
        tag_list = [tag.strip() for tag in tags.split(",")]
    
    # Get or create session
    session_id = session_manager.get_or_create_session(session_id)
    
    # Check if protein already exists
    existing = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
    if existing:
        protein_id, protein = existing
        return types.TextContent(
            type="text",
            text=f"Protein {pdb_id} already stored in session.\nSession ID: {session_id}\nProtein ID: {protein_id}\nChains: {protein.chain_count or 'N/A'}\nResidues: {protein.residue_count or 'N/A'}"
        )
    
    # Save PDB content to file if provided
    file_path = None
    if pdb_content:
        file_id = session_manager.save_file(
            session_id=session_id,
            content=pdb_content,
            filename=f"{pdb_id}.pdb",
            file_type="pdb",
            tags=tag_list
        )
        file_entry = session_manager.get_file(session_id, file_id)
        if file_entry:
            file_path = file_entry.file_path
    
    # Store protein
    protein_id = session_manager.add_protein(
        session_id=session_id,
        pdb_id=pdb_id,
        name=name,
        pdb_content=pdb_content,
        file_path=file_path,
        tags=tag_list
    )
    
    # Get stored protein for details
    protein = session_manager.get_protein(session_id, protein_id)
    
    return types.TextContent(
        type="text",
        text=f"Protein stored successfully!\nSession ID: {session_id}\nProtein ID: {protein_id}\nPDB ID: {pdb_id}\nName: {name or 'N/A'}\nChains: {protein.chain_count or 'N/A'}\nResidues: {protein.residue_count or 'N/A'}\nFile saved: {'Yes' if file_path else 'No'}\nTags: {', '.join(tag_list) if tag_list else 'None'}"
    )

def list_stored_proteins(session_id: str, tags: Optional[str] = None) -> types.TextContent:
    """List proteins stored in the session
    
    Args:
        session_id: Session ID to query
        tags: Optional comma-separated list of tags to filter by
    
    Returns:
        TextContent with list of stored proteins
    """
    
    # Parse tags if provided
    tag_filter = None
    if tags:
        tag_filter = [tag.strip() for tag in tags.split(",")]
    
    proteins = session_manager.list_proteins(session_id, tag_filter)
    
    if not proteins:
        return types.TextContent(
            type="text",
            text="No proteins found in this session."
        )
    
    # Create summary table
    lines = ["Stored Proteins:", "=" * 50]
    for prot_id, protein in proteins.items():
        lines.append(f"ID: {prot_id[:8]}...")
        lines.append(f"  PDB ID: {protein.pdb_id}")
        lines.append(f"  Name: {protein.name or 'N/A'}")
        lines.append(f"  Chains: {protein.chain_count or 'N/A'}")
        lines.append(f"  Residues: {protein.residue_count or 'N/A'}")
        lines.append(f"  File saved: {'Yes' if protein.file_path else 'No'}")
        lines.append(f"  Tags: {', '.join(protein.tags) if protein.tags else 'None'}")
        lines.append(f"  Created: {protein.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("-" * 30)
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

def convert_ligand_pdb_to_sdf(session_id: str, file_id: Optional[str] = None,
                              pdb_path: Optional[str] = None, res_name: Optional[str] = None,
                              save_to_session: bool = True) -> types.TextContent:
    """Convert a ligand PDB file to SDF format with proper bond orders and extract SMILES
    
    Args:
        session_id: Session ID
        file_id: File ID of a stored PDB file (optional)
        pdb_path: Direct path to PDB file (optional)
        res_name: Residue name of the ligand (for proper bond order assignment)
        save_to_session: Whether to save the SDF file to session
    
    Returns:
        TextContent with conversion results and SMILES
    """
    
    # Import ProDy
    try:
        import prody
    except ImportError:
        return types.TextContent(
            type="text",
            text="Error: ProDy is required for proper ligand conversion"
        )
    
    # Determine the PDB source and parse with ProDy
    original_filename = "ligand"
    ligand_atomgroup = None
    
    if file_id:
        # Get from stored file
        file_entry = session_manager.get_file(session_id, file_id)
        if file_entry:
            original_filename = file_entry.original_name.replace('.pdb', '')
            ligand_atomgroup = prody.parsePDB(file_entry.file_path)
    elif pdb_path:
        # Read from path
        try:
            original_filename = os.path.basename(pdb_path).replace('.pdb', '')
            ligand_atomgroup = prody.parsePDB(pdb_path)
        except Exception as e:
            return types.TextContent(
                type="text",
                text=f"Error reading PDB file: {str(e)}"
            )
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide either file_id or pdb_path"
        )
    
    if ligand_atomgroup is None:
        return types.TextContent(
            type="text",
            text="Error: Could not parse PDB file with ProDy"
        )
    
    # If no residue name provided, try to get it from the ligand
    if not res_name:
        res_names = set(ligand_atomgroup.getResnames())
        # Filter out common non-ligand residues
        res_names = res_names - {'HOH', 'WAT', 'H2O', 'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE'}
        if len(res_names) == 1:
            res_name = list(res_names)[0]
        elif len(res_names) > 1:
            return types.TextContent(
                type="text",
                text=f"Multiple residues found: {', '.join(res_names)}. Please specify res_name parameter."
            )
        else:
            return types.TextContent(
                type="text",
                text="No ligand residues found in PDB file"
            )
    
    results = [f"Converting {res_name} to SDF with proper bond orders", "=" * 40]
    
    # Convert to SDF with proper bond orders
    sdf_content, smiles = ligand_atomgroup_to_sdf(ligand_atomgroup, res_name)
    
    if not sdf_content or not smiles:
        # Try fallback method
        results.append("Warning: Could not convert with bond orders, trying basic conversion...")
        
        # Read PDB content for fallback
        if file_id:
            pdb_content = session_manager.read_file_content(session_id, file_id)
        else:
            with open(pdb_path, 'r') as f:
                pdb_content = f.read()
        
        sdf_content = pdb_to_sdf(pdb_content, res_name=res_name)
        
        if sdf_content:
            smiles = extract_smiles_from_sdf(sdf_content)
        
        if not sdf_content or not smiles:
            return types.TextContent(
                type="text",
                text="Error: Failed to convert PDB to SDF format"
            )
    
    results.append("Conversion successful!")
    results.append(f"SMILES: {smiles}")
    
    # Calculate properties
    properties = calculate_ligand_properties(smiles)
    if properties:
        results.append("\nMolecular Properties:")
        for prop, value in properties.items():
            if isinstance(value, float):
                results.append(f"  {prop}: {value:.3f}")
            else:
                results.append(f"  {prop}: {value}")
    
    # Save to session if requested
    if save_to_session:
        # Save SDF file
        sdf_filename = f"{original_filename}.sdf"
        sdf_file_id = session_manager.save_file(
            session_id=session_id,
            content=sdf_content,
            filename=sdf_filename,
            file_type="sdf",
            tags=["ligand", "converted", res_name]
        )
        file_entry = session_manager.get_file(session_id, sdf_file_id)
        if file_entry:
            results.append(f"\nSDF file saved: {file_entry.file_path}")
            results.append(f"File ID: {sdf_file_id}")
        
        # Store as molecule
        mol_id = session_manager.add_molecule(
            session_id=session_id,
            smiles=smiles,
            name=f"{original_filename}_{res_name}",
            properties=properties,
            tags=["from_pdb", "ligand", res_name]
        )
        results.append(f"Molecule stored with ID: {mol_id}")
    
    return types.TextContent(
        type="text",
        text="\n".join(results)
    )


def convert_protein_ligands_to_sdf(session_id: str, protein_id: Optional[str] = None,
                                   pdb_id: Optional[str] = None) -> types.TextContent:
    """Convert all ligands from a stored protein to SDF format with SMILES using proper bond orders
    
    Args:
        session_id: Session ID
        protein_id: Protein ID (optional if pdb_id provided)
        pdb_id: PDB ID (optional if protein_id provided)
    
    Returns:
        TextContent with conversion results for all ligands
    """
    
    # Get the protein entry
    if protein_id:
        protein = session_manager.get_protein(session_id, protein_id)
    elif pdb_id:
        result = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
        if result:
            protein_id, protein = result
        else:
            protein = None
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide either protein_id or pdb_id"
        )
    
    if not protein:
        return types.TextContent(
            type="text",
            text=f"Protein not found in session"
        )
    
    # Check for ligand files
    if not protein.properties or 'components' not in protein.properties:
        return types.TextContent(
            type="text",
            text=f"No extracted components found for {protein.pdb_id}"
        )
    
    components = protein.properties['components']
    
    if 'ligand_files' not in components:
        return types.TextContent(
            type="text",
            text=f"No ligand files found for {protein.pdb_id}"
        )
    
    results = [f"Converting ligands from {protein.pdb_id.upper()}", "=" * 50]
    converted_count = 0
    
    # Import ProDy for parsing ligand files
    try:
        import prody
    except ImportError:
        return types.TextContent(
            type="text",
            text="Error: ProDy is required for proper ligand conversion"
        )
    
    # Convert each ligand file
    for ligand_name, ligand_path in components['ligand_files'].items():
        results.append(f"\nProcessing {ligand_name}:")
        
        try:
            # Parse ligand PDB with ProDy
            ligand_atomgroup = prody.parsePDB(ligand_path)
            
            if ligand_atomgroup is None:
                results.append(f"  Error: Could not parse {ligand_path}")
                continue
            
            # Convert to SDF with proper bond orders
            sdf_content, smiles = ligand_atomgroup_to_sdf(ligand_atomgroup, ligand_name)
            
            if sdf_content and smiles:
                # Save SDF file
                sdf_filename = f"{protein.pdb_id}_{ligand_name}.sdf"
                sdf_file_id = session_manager.save_file(
                    session_id=session_id,
                    content=sdf_content,
                    filename=sdf_filename,
                    file_type="sdf",
                    tags=["ligand", protein.pdb_id, ligand_name]
                )
                
                file_entry = session_manager.get_file(session_id, sdf_file_id)
                if file_entry:
                    results.append(f"  SDF saved: {file_entry.file_path}")
                
                results.append(f"  SMILES: {smiles}")
                
                # Calculate properties and store as molecule
                properties = calculate_ligand_properties(smiles)
                mol_id = session_manager.add_molecule(
                    session_id=session_id,
                    smiles=smiles,
                    name=f"{protein.pdb_id}_{ligand_name}",
                    properties=properties,
                    tags=["from_pdb", protein.pdb_id, ligand_name]
                )
                results.append(f"  Stored as molecule: {mol_id[:8]}...")
                
                # Show key properties
                if properties:
                    results.append(f"  MW: {properties.get('molecular_weight', 'N/A'):.1f}, LogP: {properties.get('logp', 'N/A'):.2f}")
                
                converted_count += 1
            else:
                # Fallback to basic conversion if ProDy method fails
                results.append(f"  Warning: Could not convert with bond orders, trying basic conversion...")
                
                with open(ligand_path, 'r') as f:
                    pdb_content = f.read()
                
                sdf_content = pdb_to_sdf(pdb_content, res_name=ligand_name)
                
                if sdf_content:
                    # Save SDF file
                    sdf_filename = f"{protein.pdb_id}_{ligand_name}_basic.sdf"
                    sdf_file_id = session_manager.save_file(
                        session_id=session_id,
                        content=sdf_content,
                        filename=sdf_filename,
                        file_type="sdf",
                        tags=["ligand", protein.pdb_id, ligand_name, "basic"]
                    )
                    
                    smiles = extract_smiles_from_sdf(sdf_content)
                    if smiles:
                        results.append(f"  Basic SMILES: {smiles}")
                        converted_count += 1
                else:
                    results.append(f"  Error: Failed to convert to SDF")
                
        except Exception as e:
            results.append(f"  Error: {str(e)}")
    
    results.append(f"\n{converted_count} ligand(s) successfully converted")
    
    return types.TextContent(
        type="text",
        text="\n".join(results)
    )


def get_ligand_image(session_id: str, molecule_id: Optional[str] = None,
                    smiles: Optional[str] = None) -> Union[types.ImageContent, types.TextContent]:
    """Generate an image of a ligand molecule
    
    Args:
        session_id: Session ID
        molecule_id: ID of a stored molecule (optional if smiles provided)
        smiles: SMILES string (optional if molecule_id provided)
    
    Returns:
        ImageContent with molecule visualization or TextContent with error
    """
    from medchem_mcp_server.molecule_tools import smiles_to_image
    
    # Get SMILES
    if molecule_id:
        molecule = session_manager.get_molecule(session_id, molecule_id)
        if not molecule:
            return types.TextContent(
                type="text",
                text=f"Error: Molecule with ID {molecule_id} not found in session"
            )
        smiles = molecule.smiles
    elif not smiles:
        return types.TextContent(
            type="text",
            text="Error: Must provide either molecule_id or smiles parameter"
        )
    
    # Validate SMILES is not empty
    if not smiles or smiles.strip() == "":
        return types.TextContent(
            type="text",
            text="Error: Empty or invalid SMILES string"
        )
    
    # Generate image
    return smiles_to_image(smiles)


def prepare_protein_for_docking(session_id: str, protein_id: Optional[str] = None,
                               pdb_id: Optional[str] = None, file_id: Optional[str] = None,
                               fix_structure: bool = True, add_hydrogens: bool = True,
                               ph: float = 7.0) -> types.TextContent:
    """Prepare a protein for docking by fixing structure and converting to PDBQT format
    
    Args:
        session_id: Session ID
        protein_id: Protein ID (optional)
        pdb_id: PDB ID (optional)
        file_id: File ID of a protein PDB file (optional)
        fix_structure: Whether to fix missing atoms and residues
        add_hydrogens: Whether to add hydrogens
        ph: pH for hydrogen addition
    
    Returns:
        TextContent with preparation results
    """
    
    # Get protein file path
    protein_path = None
    protein_name = "protein"
    
    if file_id:
        file_entry = session_manager.get_file(session_id, file_id)
        if file_entry:
            protein_path = file_entry.file_path
            protein_name = file_entry.original_name.replace('.pdb', '')
    elif protein_id or pdb_id:
        # Get from stored protein
        if protein_id:
            protein = session_manager.get_protein(session_id, protein_id)
        else:
            result = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
            if result:
                protein_id, protein = result
            else:
                protein = None
        
        if protein:
            protein_name = protein.pdb_id
            # Check for extracted protein file
            if protein.properties and 'components' in protein.properties:
                if 'protein_file' in protein.properties['components']:
                    protein_path = protein.properties['components']['protein_file']
            
            # Fallback to main file
            if not protein_path and protein.file_path:
                protein_path = protein.file_path
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide protein_id, pdb_id, or file_id"
        )
    
    if not protein_path or not os.path.exists(protein_path):
        return types.TextContent(
            type="text",
            text="Error: Protein file not found"
        )
    
    results = [f"Preparing protein {protein_name} for docking", "=" * 50]
    
    try:
        # Get session temp directory
        session = session_manager.get_session(session_id)
        if not session or not session.tmp_dir:
            return types.TextContent(
                type="text",
                text="Error: Session temp directory not available"
            )
        
        # Fix structure if requested
        fixed_path = protein_path
        if fix_structure:
            results.append("Fixing protein structure...")
            fixed_filename = f"{protein_name}_fixed.pdb"
            fixed_path = os.path.join(session.tmp_dir, fixed_filename)
            
            try:
                fixed_path = add_missing_protein_atoms(
                    protein_path, 
                    fixed_path,
                    add_hydrogens=add_hydrogens,
                    ph=ph
                )
                results.append(f"  ✓ Fixed structure saved: {fixed_filename}")
                
                # Save fixed PDB to session
                with open(fixed_path, 'r') as f:
                    fixed_content = f.read()
                
                fixed_file_id = session_manager.save_file(
                    session_id=session_id,
                    content=fixed_content,
                    filename=fixed_filename,
                    file_type="pdb",
                    tags=["protein", "fixed", protein_name]
                )
                results.append(f"  ✓ Fixed PDB file ID: {fixed_file_id}")
                
            except Exception as e:
                results.append(f"  ⚠ Warning: Could not fix structure: {str(e)}")
                results.append("  Using original structure...")
        
        # Convert to PDBQT
        results.append("\nConverting to PDBQT format...")
        pdbqt_filename = f"{protein_name}.pdbqt"
        pdbqt_path = os.path.join(session.tmp_dir, pdbqt_filename)
        
        pdbqt_path = write_protein_pdbqt(fixed_path, pdbqt_path, allow_bad_res=True)
        results.append(f"  ✓ PDBQT created: {pdbqt_filename}")
        
        # Save PDBQT to session
        with open(pdbqt_path, 'r') as f:
            pdbqt_content = f.read()
        
        pdbqt_file_id = session_manager.save_file(
            session_id=session_id,
            content=pdbqt_content,
            filename=pdbqt_filename,
            file_type="pdbqt",
            tags=["protein", "docking", protein_name]
        )
        
        results.append(f"  ✓ PDBQT file ID: {pdbqt_file_id}")
        results.append(f"\n✅ Protein ready for docking!")
        results.append(f"PDBQT path: {pdbqt_path}")
        
    except Exception as e:
        results.append(f"\n❌ Error: {str(e)}")
        logger.error(f"Error preparing protein for docking: {e}", exc_info=True)
    
    return types.TextContent(
        type="text",
        text="\n".join(results)
    )


def prepare_ligand_for_docking(session_id: str, molecule_id: Optional[str] = None,
                              file_id: Optional[str] = None, smiles: Optional[str] = None,
                              name: Optional[str] = None) -> types.TextContent:
    """Prepare a ligand for docking by converting to PDBQT format
    
    Args:
        session_id: Session ID
        molecule_id: Molecule ID from session (optional)
        file_id: File ID of a ligand file (SDF/PDB) (optional)
        smiles: SMILES string (optional)
        name: Name for the ligand
    
    Returns:
        TextContent with preparation results
    """
    
    results = ["Preparing ligand for docking", "=" * 50]
    
    try:
        # Get session temp directory
        session = session_manager.get_session(session_id)
        if not session or not session.tmp_dir:
            return types.TextContent(
                type="text",
                text="Error: Session temp directory not available"
            )
        
        ligand_name = name or "ligand"
        sdf_path = None
        
        # Get or create SDF file
        if molecule_id:
            # Get molecule from session
            molecule = session_manager.get_molecule(session_id, molecule_id)
            if not molecule:
                return types.TextContent(
                    type="text",
                    text=f"Error: Molecule {molecule_id} not found"
                )
            
            ligand_name = molecule.name or ligand_name
            smiles = molecule.smiles
            results.append(f"Using molecule: {ligand_name}")
            results.append(f"SMILES: {smiles}")
            
            # Create SDF from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UFFOptimizeMolecule(mol)
                sdf_content = Chem.MolToMolBlock(mol)
                
                sdf_filename = f"{ligand_name}.sdf"
                sdf_path = os.path.join(session.tmp_dir, sdf_filename)
                with open(sdf_path, 'w') as f:
                    f.write(sdf_content)
                
        elif file_id:
            # Get file from session
            file_entry = session_manager.get_file(session_id, file_id)
            if not file_entry:
                return types.TextContent(
                    type="text",
                    text=f"Error: File {file_id} not found"
                )
            
            ligand_name = file_entry.original_name.replace('.sdf', '').replace('.pdb', '')
            
            if file_entry.file_type == 'sdf':
                sdf_path = file_entry.file_path
                results.append(f"Using SDF file: {file_entry.original_name}")
            elif file_entry.file_type == 'pdb':
                # Convert PDB to SDF first
                results.append(f"Converting PDB to SDF...")
                with open(file_entry.file_path, 'r') as f:
                    pdb_content = f.read()
                
                # Try to get residue name
                res_names = set()
                for line in pdb_content.split('\n'):
                    if line.startswith(('ATOM', 'HETATM')):
                        if len(line) > 20:
                            res_names.add(line[17:20].strip())
                
                res_names = res_names - {'HOH', 'WAT', 'H2O', 'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE'}
                res_name = list(res_names)[0] if len(res_names) == 1 else None
                
                sdf_content = pdb_to_sdf(pdb_content, res_name=res_name)
                if sdf_content:
                    sdf_filename = f"{ligand_name}.sdf"
                    sdf_path = os.path.join(session.tmp_dir, sdf_filename)
                    with open(sdf_path, 'w') as f:
                        f.write(sdf_content)
                else:
                    return types.TextContent(
                        type="text",
                        text="Error: Could not convert PDB to SDF"
                    )
                    
        elif smiles:
            # Create SDF from SMILES
            results.append(f"Creating ligand from SMILES: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UFFOptimizeMolecule(mol)
                sdf_content = Chem.MolToMolBlock(mol)
                
                sdf_filename = f"{ligand_name}.sdf"
                sdf_path = os.path.join(session.tmp_dir, sdf_filename)
                with open(sdf_path, 'w') as f:
                    f.write(sdf_content)
            else:
                return types.TextContent(
                    type="text",
                    text="Error: Invalid SMILES string"
                )
        else:
            return types.TextContent(
                type="text",
                text="Error: Must provide molecule_id, file_id, or smiles"
            )
        
        if not sdf_path or not os.path.exists(sdf_path):
            return types.TextContent(
                type="text",
                text="Error: Could not create SDF file"
            )
        
        # Convert to PDBQT
        results.append("\nConverting to PDBQT format...")
        pdbqt_filename = f"{ligand_name}.pdbqt"
        pdbqt_path = os.path.join(session.tmp_dir, pdbqt_filename)
        
        pdbqt_path = write_ligand_pdbqt(sdf_path, pdbqt_path)
        results.append(f"  ✓ PDBQT created: {pdbqt_filename}")
        
        # Save PDBQT to session
        with open(pdbqt_path, 'r') as f:
            pdbqt_content = f.read()
        
        pdbqt_file_id = session_manager.save_file(
            session_id=session_id,
            content=pdbqt_content,
            filename=pdbqt_filename,
            file_type="pdbqt",
            tags=["ligand", "docking", ligand_name]
        )
        
        results.append(f"  ✓ PDBQT file ID: {pdbqt_file_id}")
        results.append(f"\n✅ Ligand ready for docking!")
        results.append(f"PDBQT path: {pdbqt_path}")
        
    except Exception as e:
        results.append(f"\n❌ Error: {str(e)}")
        logger.error(f"Error preparing ligand for docking: {e}", exc_info=True)
    
    return types.TextContent(
        type="text",
        text="\n".join(results)
    )


def get_protein_component_files(session_id: str, protein_id: Optional[str] = None,
                               pdb_id: Optional[str] = None, component: Optional[str] = "all") -> types.TextContent:
    """Get file paths for extracted protein components
    
    Args:
        session_id: Session ID
        protein_id: Protein ID (optional if pdb_id provided)
        pdb_id: PDB ID (optional if protein_id provided)
        component: Component to get ('protein', 'ligands', 'waters', 'ions', or 'all')
    
    Returns:
        TextContent with file paths for requested components
    """
    
    # Get the protein entry
    if protein_id:
        protein = session_manager.get_protein(session_id, protein_id)
    elif pdb_id:
        result = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
        if result:
            protein_id, protein = result
        else:
            protein = None
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide either protein_id or pdb_id"
        )
    
    if not protein:
        return types.TextContent(
            type="text",
            text=f"Protein not found in session"
        )
    
    if not protein.properties or 'components' not in protein.properties:
        return types.TextContent(
            type="text",
            text=f"No extracted components found for {protein.pdb_id}. The protein may not have been processed with component extraction."
        )
    
    components = protein.properties['components']
    lines = [f"Component files for {protein.pdb_id.upper()}:", "=" * 50]
    
    if component == "all" or component == "protein":
        if 'protein_file' in components:
            lines.append(f"Protein only: {components['protein_file']}")
            lines.append(f"  Atoms: {components.get('protein_atoms', 'N/A')}")
    
    if component == "all" or component == "ligands":
        if 'all_ligands_file' in components:
            lines.append(f"All ligands: {components['all_ligands_file']}")
            lines.append(f"  Atoms: {components.get('ligand_atoms', 'N/A')}")
            lines.append(f"  Types: {', '.join(components.get('ligand_names', []))}")
        
        if 'ligand_files' in components:
            lines.append("Individual ligand files:")
            for lig_name, lig_path in components['ligand_files'].items():
                lines.append(f"  {lig_name}: {lig_path}")
    
    if component == "all" or component == "waters":
        if 'waters_file' in components:
            lines.append(f"Waters: {components['waters_file']}")
            lines.append(f"  Molecules: {components.get('water_molecules', 'N/A')}")
    
    if component == "all" or component == "ions":
        if 'ions_file' in components:
            lines.append(f"Ions: {components['ions_file']}")
            lines.append(f"  Atoms: {components.get('ion_atoms', 'N/A')}")
            lines.append(f"  Types: {', '.join(components.get('ion_types', []))}")
    
    if len(lines) == 2:  # Only header lines
        lines.append(f"No {component} component files found.")
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

def analyze_stored_protein(session_id: str, protein_id: Optional[str] = None,
                          pdb_id: Optional[str] = None) -> types.TextContent:
    """Analyze a stored protein using ProDy utilities to segment components
    
    Args:
        session_id: Session ID
        protein_id: Protein ID to analyze (optional if pdb_id provided)
        pdb_id: PDB ID to analyze (optional if protein_id provided)
    
    Returns:
        TextContent with detailed protein analysis
    """
    
    # Get the protein entry
    if protein_id:
        protein = session_manager.get_protein(session_id, protein_id)
    elif pdb_id:
        result = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
        if result:
            protein_id, protein = result
        else:
            protein = None
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide either protein_id or pdb_id"
        )
    
    if not protein or not protein.file_path:
        return types.TextContent(
            type="text",
            text=f"Protein not found in session or has no saved file"
        )
    
    if not UTILITIES_AVAILABLE:
        return types.TextContent(
            type="text",
            text="Advanced protein analysis not available. Please install scipy and ensure all dependencies are met."
        )
    
    try:
        import prody
        import os
        
        # Parse the PDB file using ProDy
        if not os.path.exists(protein.file_path):
            return types.TextContent(
                type="text",
                text=f"PDB file not found at {protein.file_path}"
            )
        
        structure = prody.parsePDB(protein.file_path)
        
        if structure is None:
            return types.TextContent(
                type="text",
                text=f"Failed to parse PDB file for {protein.pdb_id}"
            )
        
        # Analyze components using utilities
        protein_atoms = get_protein(structure)
        ligands = get_ligands(structure)
        waters = get_waters(structure)
        ions = get_ions(structure)
        nucleic_acids = get_nucleic_acids(structure)
        non_standard_aa = get_non_standard_amino_acids(structure)
        
        # Build analysis report
        lines = [
            f"Protein Structure Analysis for {protein.pdb_id.upper()}",
            "=" * 50,
            f"Total atoms: {structure.numAtoms()}",
            ""
        ]
        
        # Protein component
        if protein_atoms:
            lines.append(f"Protein Component:")
            lines.append(f"  Atoms: {len(protein_atoms)}")
            chains = set(protein_atoms.getChids())
            lines.append(f"  Chains: {', '.join(sorted(chains))}")
            
            # Count residues per chain
            for chain in sorted(chains):
                chain_sel = protein_atoms.select(f'chain {chain}').toAtomGroup()
                if chain_sel:
                    chain_residues = set()
                    for res in chain_sel.iterResidues():
                        chain_residues.add(res.getResnum())
                    lines.append(f"    Chain {chain}: {len(chain_residues)} residues")
            lines.append("")
        
        # Ligands
        if ligands:
            lines.append(f"Ligands:")
            lines.append(f"  Total atoms: {len(ligands)}")
            ligand_names = set(ligands.getResnames())
            lines.append(f"  Molecules: {', '.join(sorted(ligand_names))}")
            
            # Count each ligand type
            for lig_name in sorted(ligand_names):
                lig_sel = ligands.select(f'resname {lig_name}').toAtomGroup()
                if lig_sel:
                    # Count unique residue numbers for this ligand
                    unique_ligs = set()
                    for res in lig_sel.iterResidues():
                        unique_ligs.add((res.getChid(), res.getResnum()))
                    lines.append(f"    {lig_name}: {len(unique_ligs)} instance(s)")
            lines.append("")
        
        # Waters
        if waters:
            lines.append(f"Water Molecules:")
            lines.append(f"  Total atoms: {len(waters)}")
            water_residues = set()
            for res in waters.iterResidues():
                water_residues.add((res.getChid(), res.getResnum()))
            lines.append(f"  Water molecules: {len(water_residues)}")
            lines.append("")
        
        # Ions
        if ions:
            lines.append(f"Ions:")
            lines.append(f"  Total atoms: {len(ions)}")
            ion_names = set(ions.getResnames())
            lines.append(f"  Ion types: {', '.join(sorted(ion_names))}")
            
            # Count each ion type
            for ion_name in sorted(ion_names):
                ion_sel = ions.select(f'resname {ion_name}').toAtomGroup()
                if ion_sel:
                    lines.append(f"    {ion_name}: {len(ion_sel)} atom(s)")
            lines.append("")
        
        # Nucleic acids
        if nucleic_acids:
            lines.append(f"Nucleic Acids:")
            lines.append(f"  Total atoms: {len(nucleic_acids)}")
            na_chains = set(nucleic_acids.getChids())
            lines.append(f"  Chains: {', '.join(sorted(na_chains))}")
            lines.append("")
        
        # Non-standard amino acids
        if non_standard_aa:
            lines.append(f"Non-standard Amino Acids:")
            lines.append(f"  Total atoms: {len(non_standard_aa)}")
            ns_names = set(non_standard_aa.getResnames())
            lines.append(f"  Types: {', '.join(sorted(ns_names))}")
            
            # Count each type
            for ns_name in sorted(ns_names):
                ns_sel = non_standard_aa.select(f'resname {ns_name}').toAtomGroup()
                if ns_sel:
                    unique_ns = set()
                    for res in ns_sel.iterResidues():
                        unique_ns.add((res.getChid(), res.getResnum()))
                    lines.append(f"    {ns_name}: {len(unique_ns)} residue(s)")
            lines.append("")
        
        # Summary
        lines.append("Summary:")
        lines.append(f"  Protein atoms: {len(protein_atoms) if protein_atoms else 0}")
        lines.append(f"  Ligand atoms: {len(ligands) if ligands else 0}")
        lines.append(f"  Water atoms: {len(waters) if waters else 0}")
        lines.append(f"  Ion atoms: {len(ions) if ions else 0}")
        lines.append(f"  Nucleic acid atoms: {len(nucleic_acids) if nucleic_acids else 0}")
        lines.append(f"  Non-standard AA atoms: {len(non_standard_aa) if non_standard_aa else 0}")
        
        return types.TextContent(
            type="text",
            text="\n".join(lines)
        )
        
    except Exception as e:
        logger.error(f"Error analyzing protein: {e}", exc_info=True)
        return types.TextContent(
            type="text",
            text=f"Error analyzing protein: {str(e)}"
        )

def get_stored_protein(session_id: str, protein_id: Optional[str] = None, 
                      pdb_id: Optional[str] = None) -> types.TextContent:
    """Get details of a specific stored protein
    
    Args:
        session_id: Session ID
        protein_id: Protein ID to retrieve (optional if pdb_id provided)
        pdb_id: PDB ID to retrieve (optional if protein_id provided)
    
    Returns:
        TextContent with protein details
    """
    
    if protein_id:
        protein = session_manager.get_protein(session_id, protein_id)
    elif pdb_id:
        result = session_manager.get_protein_by_pdb_id(session_id, pdb_id)
        if result:
            protein_id, protein = result
        else:
            protein = None
    else:
        return types.TextContent(
            type="text",
            text="Error: Must provide either protein_id or pdb_id"
        )
    
    if not protein:
        return types.TextContent(
            type="text",
            text=f"Protein not found in session {session_id}"
        )
    
    details = [
        f"Protein Details:",
        f"ID: {protein_id}",
        f"PDB ID: {protein.pdb_id}",
        f"Name: {protein.name or 'N/A'}",
        f"Chains: {protein.chain_count or 'N/A'}",
        f"Residues: {protein.residue_count or 'N/A'}",
        f"File path: {protein.file_path or 'Not saved'}",
        f"Tags: {', '.join(protein.tags) if protein.tags else 'None'}",
        f"Created: {protein.created_at.strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    
    if protein.properties:
        details.append("")
        details.append("Properties:")
        for prop, value in protein.properties.items():
            if prop == 'chain_info' and isinstance(value, dict):
                details.append(f"  Chain Information:")
                for chain_id, res_count in value.items():
                    details.append(f"    Chain {chain_id}: {res_count} residues")
            elif prop == 'components' and isinstance(value, dict):
                details.append(f"  Components:")
                for comp_name, comp_value in value.items():
                    if isinstance(comp_value, list):
                        details.append(f"    {comp_name}: {', '.join(map(str, comp_value))}")
                    elif isinstance(comp_value, dict):
                        # Handle nested dicts like ligand_files
                        details.append(f"    {comp_name}:")
                        for sub_key, sub_value in comp_value.items():
                            details.append(f"      {sub_key}: {sub_value}")
                    elif comp_name.endswith('_file') or comp_name.endswith('_files'):
                        # Show file paths
                        details.append(f"    {comp_name}: {comp_value}")
                    else:
                        details.append(f"    {comp_name}: {comp_value}")
            else:
                details.append(f"  {prop}: {value}")
    
    return types.TextContent(
        type="text",
        text="\n".join(details)
    )

# File management tools

def list_stored_files(session_id: str, file_type: Optional[str] = None, 
                     tags: Optional[str] = None) -> types.TextContent:
    """List files stored in the session
    
    Args:
        session_id: Session ID to query
        file_type: Optional file type to filter by (e.g., 'pdb', 'sdf')
        tags: Optional comma-separated list of tags to filter by
    
    Returns:
        TextContent with list of stored files
    """
    
    # Parse tags if provided
    tag_filter = None
    if tags:
        tag_filter = [tag.strip() for tag in tags.split(",")]
    
    files = session_manager.list_files(session_id, file_type, tag_filter)
    
    if not files:
        return types.TextContent(
            type="text",
            text="No files found in this session."
        )
    
    # Create summary table
    lines = ["Stored Files:", "=" * 50]
    for file_id, file_entry in files.items():
        lines.append(f"ID: {file_id[:8]}...")
        lines.append(f"  Name: {file_entry.original_name}")
        lines.append(f"  Type: {file_entry.file_type}")
        lines.append(f"  Size: {file_entry.size_bytes} bytes ({file_entry.size_bytes / 1024:.2f} KB)")
        lines.append(f"  Path: {file_entry.file_path}")
        lines.append(f"  Tags: {', '.join(file_entry.tags) if file_entry.tags else 'None'}")
        lines.append(f"  Created: {file_entry.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("-" * 30)
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

def get_stored_file(session_id: str, file_id: str) -> types.TextContent:
    """Get details of a specific stored file
    
    Args:
        session_id: Session ID
        file_id: File ID to retrieve
    
    Returns:
        TextContent with file details
    """
    
    file_entry = session_manager.get_file(session_id, file_id)
    
    if not file_entry:
        return types.TextContent(
            type="text",
            text=f"File with ID {file_id} not found in session {session_id}"
        )
    
    details = [
        f"File Details:",
        f"ID: {file_id}",
        f"Original name: {file_entry.original_name}",
        f"Type: {file_entry.file_type}",
        f"Size: {file_entry.size_bytes} bytes ({file_entry.size_bytes / 1024:.2f} KB)",
        f"Path: {file_entry.file_path}",
        f"Tags: {', '.join(file_entry.tags) if file_entry.tags else 'None'}",
        f"Created: {file_entry.created_at.strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    
    if file_entry.properties:
        details.append("")
        details.append("Properties:")
        for prop, value in file_entry.properties.items():
            details.append(f"  {prop}: {value}")
    
    return types.TextContent(
        type="text",
        text="\n".join(details)
    )

def read_stored_file(session_id: str, file_id: str, max_lines: Optional[int] = 100) -> types.TextContent:
    """Read content of a stored file
    
    Args:
        session_id: Session ID
        file_id: File ID to read
        max_lines: Maximum number of lines to return (default 100)
    
    Returns:
        TextContent with file content
    """
    
    content = session_manager.read_file_content(session_id, file_id)
    
    if content is None:
        return types.TextContent(
            type="text",
            text=f"File with ID {file_id} not found or could not be read"
        )
    
    lines = content.split('\n')
    if max_lines and len(lines) > max_lines:
        truncated_content = '\n'.join(lines[:max_lines])
        return types.TextContent(
            type="text",
            text=f"File content (first {max_lines} lines of {len(lines)} total):\n\n{truncated_content}\n\n[Truncated...]"
        )
    
    return types.TextContent(
        type="text",
        text=f"File content:\n\n{content}"
    )

# Collection management tools

def create_collection(session_id: str, collection_name: str, 
                     entity_ids: str) -> types.TextContent:
    """Create a named collection of entities (molecules, proteins, or files)
    
    Args:
        session_id: Session ID
        collection_name: Name for the collection
        entity_ids: Comma-separated list of entity IDs
    
    Returns:
        TextContent with collection creation confirmation
    """
    
    # Parse entity IDs
    ids = [eid.strip() for eid in entity_ids.split(",")]
    
    success = session_manager.create_collection(session_id, collection_name, ids)
    
    if success:
        return types.TextContent(
            type="text",
            text=f"Collection '{collection_name}' created with {len(ids)} entities"
        )
    else:
        return types.TextContent(
            type="text",
            text=f"Failed to create collection. Session {session_id} not found."
        )

def get_collection(session_id: str, collection_name: str) -> types.TextContent:
    """Get entities from a named collection
    
    Args:
        session_id: Session ID
        collection_name: Name of the collection
    
    Returns:
        TextContent with collection contents
    """
    
    entities = session_manager.get_collection(session_id, collection_name)
    
    if not entities:
        return types.TextContent(
            type="text",
            text=f"Collection '{collection_name}' not found or is empty"
        )
    
    lines = [f"Collection '{collection_name}':", "=" * 50]
    
    molecules = []
    proteins = []
    files = []
    
    for entity_id, entity_info in entities.items():
        if entity_info['type'] == 'molecule':
            molecules.append((entity_id, entity_info['data']))
        elif entity_info['type'] == 'protein':
            proteins.append((entity_id, entity_info['data']))
        elif entity_info['type'] == 'file':
            files.append((entity_id, entity_info['data']))
    
    if molecules:
        lines.append("\nMolecules:")
        for mol_id, mol in molecules:
            lines.append(f"  - {mol_id[:8]}... : {mol.name or mol.smiles[:30]}")
    
    if proteins:
        lines.append("\nProteins:")
        for prot_id, prot in proteins:
            lines.append(f"  - {prot_id[:8]}... : {prot.pdb_id} ({prot.name or 'N/A'})")
    
    if files:
        lines.append("\nFiles:")
        for file_id, file_entry in files:
            lines.append(f"  - {file_id[:8]}... : {file_entry.original_name}")
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

# Session management tools

def get_session_summary(session_id: str) -> types.TextContent:
    """Get summary of the current session state
    
    Args:
        session_id: Session ID to summarize
    
    Returns:
        TextContent with session summary
    """
    
    summary = session_manager.get_session_summary(session_id)
    
    if not summary:
        return types.TextContent(
            type="text",
            text=f"Session {session_id} not found"
        )
    
    lines = [
        f"Session Summary:",
        f"Session ID: {summary['session_id']}",
        f"Molecules stored: {summary['molecule_count']}",
        f"Proteins stored: {summary['protein_count']}",
        f"Files stored: {summary['file_count']}",
        f"Total file size: {summary['total_file_size_mb']} MB",
        f"Collections: {summary['collection_count']}",
        f"Collection names: {', '.join(summary['collections']) if summary['collections'] else 'None'}",
        f"Temp directory: {summary['tmp_directory']}",
        f"Created: {summary['created_at']}",
        f"Last accessed: {summary['last_accessed']}"
    ]
    
    return types.TextContent(
        type="text",
        text="\n".join(lines)
    )

def clear_session(session_id: str) -> types.TextContent:
    """Clear all data from the session
    
    Args:
        session_id: Session ID to clear
    
    Returns:
        TextContent with clear confirmation
    """
    
    success = session_manager.clear_session(session_id)
    
    if success:
        return types.TextContent(
            type="text",
            text=f"Session {session_id} cleared successfully. All temporary files deleted."
        )
    else:
        return types.TextContent(
            type="text",
            text=f"Session {session_id} not found"
        )

def create_session() -> types.TextContent:
    """Create a new session and return its ID
    
    Returns:
        TextContent with new session ID
    """
    
    session_id = session_manager.get_or_create_session()
    
    return types.TextContent(
        type="text",
        text=f"New session created!\nSession ID: {session_id}\n\nUse this ID for all subsequent operations in this session."
    )
