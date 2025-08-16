"""
Integration tools for ZINC fetching and molecule editing.

This module provides high-level tools that combine ZINC database access
with molecule editing operations for medicinal chemistry workflows.
"""

import logging
import json
from typing import Optional, Dict, List, Any
import mcp.types as types
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw

# Import ZINC tools
from medchem_mcp_server.zinc_tooling import (
    fetch_zinc_molecule,
    prepare_zinc_molecule_for_editing,
    standardize_molecule,
    calculate_properties
)

# Import editing tools
from medchem_mcp_server.molecule_editing import (
    load_and_prepare_molecule,
    get_bemis_murcko_scaffold,
    label_r_groups,
    create_templated_molecule,
    replace_rgroup,
    substructure_replace,
    scaffold_swap,
    apply_edit_plan,
    validate_molecule,
    fragment_molecule,
    enumerate_r_groups,
    COMMON_R_GROUPS,
    COMMON_SCAFFOLDS,
    COMMON_TRANSFORMATIONS
)

# Import session management
from medchem_mcp_server.session_state import session_manager
from medchem_mcp_server.utilities import pil_image_to_base64

logger = logging.getLogger("medchem-mcp-server")


# ============================================================================
# Main Integration Tools for MCP
# ============================================================================

def fetch_and_edit_zinc_molecule(zinc_id: str, edits_json: str, 
                                session_id: Optional[str] = None) -> types.TextContent:
    """
    Fetch a molecule from ZINC and apply edit operations.
    
    Args:
        zinc_id: ZINC identifier
        edits_json: JSON string with list of edit operations
        session_id: Optional session ID for storing results
    
    Returns:
        TextContent with edited molecule and validation results
    """
    try:
        # Fetch molecule from ZINC
        fetch_result = fetch_zinc_molecule(zinc_id)
        fetch_data = json.loads(fetch_result.text)
        
        if "error" in fetch_data or "Error" in fetch_result.text:
            return fetch_result
        
        original_smiles = fetch_data["smiles"]
        
        # Load molecule
        mol = load_and_prepare_molecule(original_smiles)
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Could not process molecule {zinc_id}"
            )
        
        # Apply edits
        edited_mol, messages = apply_edit_plan(mol, edits_json)
        
        if edited_mol is None:
            return types.TextContent(
                type="text",
                text=json.dumps({
                    "zinc_id": zinc_id,
                    "original_smiles": original_smiles,
                    "error": "Edit plan failed",
                    "messages": messages
                }, indent=2)
            )
        
        # Get edited SMILES
        edited_smiles = Chem.MolToSmiles(edited_mol)
        
        # Calculate properties for edited molecule
        edited_properties = calculate_properties(edited_mol)
        
        # Validate edited molecule
        validation = validate_molecule(edited_mol)
        
        # Store in session if requested
        stored_id = None
        if session_id:
            session_id = session_manager.get_or_create_session(session_id)
            stored_id = session_manager.add_molecule(
                session_id=session_id,
                smiles=edited_smiles,
                name=f"Edited_{zinc_id}",
                properties=edited_properties,
                tags=["zinc", "edited", zinc_id]
            )
        
        result = {
            "zinc_id": zinc_id,
            "original_smiles": original_smiles,
            "edited_smiles": edited_smiles,
            "edit_messages": messages,
            "original_properties": fetch_data.get("properties", {}),
            "edited_properties": edited_properties,
            "validation": validation,
            "stored_molecule_id": stored_id,
            "session_id": session_id
        }
        
        return types.TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error in fetch_and_edit_zinc_molecule: {e}")
        return types.TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )


def create_zinc_analog_series(zinc_id: str, r_group_type: str = "aromatic",
                             r_position: int = 1, 
                             session_id: Optional[str] = None) -> types.TextContent:
    """
    Create an analog series by systematically replacing R-groups.
    
    Args:
        zinc_id: ZINC identifier of parent molecule
        r_group_type: Type of R-groups to enumerate (alkyl, aromatic, etc.)
        r_position: Which R-group position to modify
        session_id: Optional session ID for storing results
    
    Returns:
        TextContent with analog series information
    """
    try:
        # Fetch parent molecule
        fetch_result = fetch_zinc_molecule(zinc_id)
        fetch_data = json.loads(fetch_result.text)
        
        if "error" in fetch_data:
            return fetch_result
        
        parent_smiles = fetch_data["smiles"]
        parent_mol = load_and_prepare_molecule(parent_smiles)
        
        if parent_mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Could not process molecule {zinc_id}"
            )
        
        # Create template with labeled R-groups
        template_info = create_templated_molecule(parent_mol)
        
        if "error" in template_info:
            return types.TextContent(
                type="text",
                text=json.dumps(template_info, indent=2)
            )
        
        # Get R-group library
        r_groups = COMMON_R_GROUPS.get(r_group_type, COMMON_R_GROUPS["aromatic"])
        
        # Generate analogs
        analogs = []
        scaffold = get_bemis_murcko_scaffold(parent_mol)
        
        for i, fragment in enumerate(r_groups):
            try:
                # Create molecule with scaffold and fragment
                # First, need to prepare a molecule with the scaffold and dummy atoms
                scaffold_with_dummy = Chem.MolFromSmiles(f"{Chem.MolToSmiles(scaffold)}[*:{r_position}]")
                if scaffold_with_dummy is None:
                    continue
                
                # Replace the dummy with the fragment
                analog_mol = replace_rgroup(scaffold_with_dummy, r_position, fragment)
                
                if analog_mol:
                    analog_smiles = Chem.MolToSmiles(analog_mol)
                    props = calculate_properties(analog_mol)
                    validation = validate_molecule(analog_mol)
                    
                    analog_info = {
                        "index": i + 1,
                        "smiles": analog_smiles,
                        "fragment": fragment,
                        "properties": props,
                        "valid": validation["valid"],
                        "warnings": validation.get("warnings", [])
                    }
                    
                    # Store if session provided and molecule is valid
                    if session_id and validation["valid"]:
                        session_id = session_manager.get_or_create_session(session_id)
                        mol_id = session_manager.add_molecule(
                            session_id=session_id,
                            smiles=analog_smiles,
                            name=f"Analog_{zinc_id}_{i+1}",
                            properties=props,
                            tags=["zinc", "analog", zinc_id, r_group_type]
                        )
                        analog_info["stored_id"] = mol_id
                    
                    analogs.append(analog_info)
                    
            except Exception as e:
                logger.warning(f"Failed to create analog with fragment {fragment}: {e}")
                continue
        
        result = {
            "parent_zinc_id": zinc_id,
            "parent_smiles": parent_smiles,
            "scaffold_smiles": Chem.MolToSmiles(scaffold) if scaffold else None,
            "r_group_type": r_group_type,
            "r_position": r_position,
            "num_analogs": len(analogs),
            "analogs": analogs,
            "session_id": session_id
        }
        
        return types.TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error creating analog series: {e}")
        return types.TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )


def scaffold_hop_zinc_molecule(zinc_id: str, new_scaffold_type: str = "pyridyl",
                              session_id: Optional[str] = None) -> types.TextContent:
    """
    Perform scaffold hopping on a ZINC molecule.
    
    Args:
        zinc_id: ZINC identifier
        new_scaffold_type: Type of scaffold to swap to
        session_id: Optional session ID
    
    Returns:
        TextContent with scaffold-hopped molecule
    """
    try:
        # Fetch molecule
        fetch_result = fetch_zinc_molecule(zinc_id)
        fetch_data = json.loads(fetch_result.text)
        
        if "error" in fetch_data:
            return fetch_result
        
        original_smiles = fetch_data["smiles"]
        mol = load_and_prepare_molecule(original_smiles)
        
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Could not process molecule {zinc_id}"
            )
        
        # Get current scaffold
        scaffold = get_bemis_murcko_scaffold(mol)
        if scaffold is None:
            return types.TextContent(
                type="text",
                text="Error: Could not extract scaffold"
            )
        
        # Simplified scaffold hopping - replace phenyl with selected scaffold
        # This is a basic example - real scaffold hopping would be more sophisticated
        
        old_scaffold_smarts = "[*:1]c1ccc([*:2])cc1"  # phenyl pattern
        new_scaffold_smiles = COMMON_SCAFFOLDS.get(new_scaffold_type, COMMON_SCAFFOLDS["pyridyl"])
        
        # Attempt scaffold swap
        hopped_mol = scaffold_swap(mol, old_scaffold_smarts, new_scaffold_smiles)
        
        if hopped_mol is None:
            # Try simpler substructure replacement
            results = substructure_replace(mol, "c1ccccc1", "c1ncccn1", False)
            if results:
                hopped_mol = results[0]
            else:
                return types.TextContent(
                    type="text",
                    text="Error: Could not perform scaffold hopping - no matching scaffold pattern"
                )
        
        hopped_smiles = Chem.MolToSmiles(hopped_mol)
        hopped_properties = calculate_properties(hopped_mol)
        validation = validate_molecule(hopped_mol)
        
        # Store if session provided
        stored_id = None
        if session_id:
            session_id = session_manager.get_or_create_session(session_id)
            stored_id = session_manager.add_molecule(
                session_id=session_id,
                smiles=hopped_smiles,
                name=f"ScaffoldHop_{zinc_id}_{new_scaffold_type}",
                properties=hopped_properties,
                tags=["zinc", "scaffold_hop", zinc_id, new_scaffold_type]
            )
        
        result = {
            "zinc_id": zinc_id,
            "original_smiles": original_smiles,
            "original_scaffold": Chem.MolToSmiles(scaffold),
            "new_scaffold_type": new_scaffold_type,
            "hopped_smiles": hopped_smiles,
            "original_properties": fetch_data.get("properties", {}),
            "hopped_properties": hopped_properties,
            "validation": validation,
            "stored_molecule_id": stored_id,
            "session_id": session_id
        }
        
        return types.TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error in scaffold hopping: {e}")
        return types.TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )


def fragment_and_recombine_zinc(zinc_ids: str, method: str = "brics",
                               session_id: Optional[str] = None) -> types.TextContent:
    """
    Fragment multiple ZINC molecules and suggest recombinations.
    
    Args:
        zinc_ids: Comma-separated ZINC IDs
        method: Fragmentation method (brics or recap)
        session_id: Optional session ID
    
    Returns:
        TextContent with fragments and suggested recombinations
    """
    try:
        ids = [id.strip() for id in zinc_ids.split(",")]
        all_fragments = set()
        molecules = []
        
        # Collect fragments from all molecules
        for zinc_id in ids:
            fetch_result = fetch_zinc_molecule(zinc_id)
            fetch_data = json.loads(fetch_result.text)
            
            if "error" in fetch_data:
                continue
            
            smiles = fetch_data["smiles"]
            mol = load_and_prepare_molecule(smiles)
            
            if mol:
                molecules.append({
                    "zinc_id": zinc_id,
                    "smiles": smiles,
                    "mol": mol
                })
                
                fragments = fragment_molecule(mol, method)
                all_fragments.update(fragments)
        
        # Generate some recombinations (simplified example)
        recombinations = []
        fragment_list = list(all_fragments)[:10]  # Limit to first 10 fragments
        
        # Try combining pairs of fragments
        for i in range(min(5, len(fragment_list))):
            for j in range(i+1, min(len(fragment_list), i+3)):
                try:
                    # Simple concatenation - in practice you'd use reaction rules
                    frag1 = fragment_list[i]
                    frag2 = fragment_list[j]
                    
                    # Remove dummy atoms and connect
                    # This is simplified - real implementation would use BRICS.BRICSBuild
                    combined_smiles = frag1.replace("[*]", "") + frag2.replace("[*]", "")
                    combined_mol = Chem.MolFromSmiles(combined_smiles)
                    
                    if combined_mol:
                        props = calculate_properties(combined_mol)
                        validation = validate_molecule(combined_mol)
                        
                        if validation["valid"]:
                            recombinations.append({
                                "smiles": Chem.MolToSmiles(combined_mol),
                                "fragments_used": [frag1, frag2],
                                "properties": props
                            })
                            
                except Exception:
                    continue
        
        result = {
            "input_molecules": [{"zinc_id": m["zinc_id"], "smiles": m["smiles"]} 
                              for m in molecules],
            "fragmentation_method": method,
            "total_fragments": len(all_fragments),
            "unique_fragments": list(all_fragments)[:20],  # Show first 20
            "suggested_recombinations": recombinations[:10],  # Show first 10
            "session_id": session_id
        }
        
        # Store recombinations if session provided
        if session_id and recombinations:
            session_id = session_manager.get_or_create_session(session_id)
            for i, recomb in enumerate(recombinations[:5]):
                session_manager.add_molecule(
                    session_id=session_id,
                    smiles=recomb["smiles"],
                    name=f"Recombination_{i+1}",
                    properties=recomb["properties"],
                    tags=["zinc", "recombination", method]
                )
        
        return types.TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error in fragment and recombine: {e}")
        return types.TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )


def apply_medchem_transforms(zinc_id: str, transform_names: str,
                            session_id: Optional[str] = None) -> types.TextContent:
    """
    Apply common medicinal chemistry transformations to a ZINC molecule.
    
    Args:
        zinc_id: ZINC identifier
        transform_names: Comma-separated transformation names
        session_id: Optional session ID
    
    Returns:
        TextContent with transformed molecules
    """
    try:
        # Fetch molecule
        fetch_result = fetch_zinc_molecule(zinc_id)
        fetch_data = json.loads(fetch_result.text)
        
        if "error" in fetch_data:
            return fetch_result
        
        original_smiles = fetch_data["smiles"]
        mol = load_and_prepare_molecule(original_smiles)
        
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Could not process molecule {zinc_id}"
            )
        
        # Parse transformation names
        transforms = [t.strip() for t in transform_names.split(",")]
        
        results = []
        for transform_name in transforms:
            if transform_name not in COMMON_TRANSFORMATIONS:
                results.append({
                    "transform": transform_name,
                    "error": "Unknown transformation"
                })
                continue
            
            smirks = COMMON_TRANSFORMATIONS[transform_name]
            
            try:
                # Apply transformation
                from medchem_mcp_server.molecule_editing import apply_transformation
                transformed_mol = apply_transformation(mol, smirks)
                
                if transformed_mol:
                    transformed_smiles = Chem.MolToSmiles(transformed_mol)
                    props = calculate_properties(transformed_mol)
                    validation = validate_molecule(transformed_mol)
                    
                    result_entry = {
                        "transform": transform_name,
                        "smirks": smirks,
                        "smiles": transformed_smiles,
                        "properties": props,
                        "valid": validation["valid"],
                        "warnings": validation.get("warnings", [])
                    }
                    
                    # Store if valid and session provided
                    if session_id and validation["valid"]:
                        session_id = session_manager.get_or_create_session(session_id)
                        mol_id = session_manager.add_molecule(
                            session_id=session_id,
                            smiles=transformed_smiles,
                            name=f"Transform_{zinc_id}_{transform_name}",
                            properties=props,
                            tags=["zinc", "transform", zinc_id, transform_name]
                        )
                        result_entry["stored_id"] = mol_id
                    
                    results.append(result_entry)
                else:
                    results.append({
                        "transform": transform_name,
                        "error": "Transformation did not match molecule"
                    })
                    
            except Exception as e:
                results.append({
                    "transform": transform_name,
                    "error": str(e)
                })
        
        output = {
            "zinc_id": zinc_id,
            "original_smiles": original_smiles,
            "original_properties": fetch_data.get("properties", {}),
            "transformations": results,
            "session_id": session_id
        }
        
        return types.TextContent(
            type="text",
            text=json.dumps(output, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error applying transformations: {e}")
        return types.TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )


def visualize_zinc_edits(zinc_id: str, edits_json: str) -> types.ImageContent:
    """
    Fetch a ZINC molecule, apply edits, and return visualization.
    
    Args:
        zinc_id: ZINC identifier
        edits_json: JSON string with edit operations
    
    Returns:
        ImageContent with before/after visualization
    """
    try:
        # Fetch molecule
        fetch_result = fetch_zinc_molecule(zinc_id)
        fetch_data = json.loads(fetch_result.text)
        
        if "error" in fetch_data:
            # Return error as text since we expect an image
            mol = Chem.MolFromSmiles("C")  # Dummy molecule
            img = Draw.MolToImage(mol)
            img_str = pil_image_to_base64(img)
            return types.ImageContent(
                type="image", 
                data=img_str, 
                mimeType="image/png"
            )
        
        original_smiles = fetch_data["smiles"]
        original_mol = load_and_prepare_molecule(original_smiles)
        
        # Apply edits
        edited_mol, _ = apply_edit_plan(original_mol, edits_json)
        
        if edited_mol is None:
            edited_mol = original_mol  # Show original if edit failed
        
        # Create side-by-side visualization
        mols = [original_mol, edited_mol]
        legends = ["Original", "Edited"]
        
        img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(300, 300),
                                  legends=legends)
        
        img_str = pil_image_to_base64(img)
        return types.ImageContent(
            type="image", 
            data=img_str, 
            mimeType="image/png"
        )
        
    except Exception as e:
        logger.error(f"Error visualizing edits: {e}")
        # Return dummy image on error
        mol = Chem.MolFromSmiles("C")
        img = Draw.MolToImage(mol)
        img_str = pil_image_to_base64(img)
        return types.ImageContent(
            type="image", 
            data=img_str, 
            mimeType="image/png"
        )
