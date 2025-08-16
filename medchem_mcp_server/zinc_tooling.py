"""
ZINC database fetching and integration tools.

This module provides functionality to fetch molecules from the ZINC database
and prepare them for editing operations.
"""

import logging
import requests
from typing import Optional, Dict, List, Any
import mcp.types as types
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import json
import io

logger = logging.getLogger("medchem-mcp-server")

# ZINC API endpoints
ZINC_BASE_URL = "https://zinc.docking.org"
ZINC_API_BASE = f"{ZINC_BASE_URL}/api/v1"


def fetch_zinc_molecule(zinc_id: str, add_hydrogens: bool = False) -> types.TextContent:
    """
    Fetch a molecule from ZINC database by ZINC ID.
    
    Args:
        zinc_id: ZINC identifier (e.g., "ZINC000000000001")
        add_hydrogens: Whether to add explicit hydrogens
    
    Returns:
        TextContent with molecule information including SMILES, properties, and status
    """
    try:
        # Clean up ZINC ID (remove ZINC prefix if provided with lowercase)
        zinc_id = zinc_id.upper()
        if not zinc_id.startswith("ZINC"):
            zinc_id = f"ZINC{zinc_id}"
        
        # Fetch molecule data from ZINC
        url = f"{ZINC_API_BASE}/substances/{zinc_id}"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 404:
            return types.TextContent(
                type="text",
                text=f"Error: Molecule {zinc_id} not found in ZINC database"
            )
        
        response.raise_for_status()
        data = response.json()
        
        # Extract SMILES
        smiles = data.get("smiles", None)
        if not smiles:
            return types.TextContent(
                type="text",
                text=f"Error: No SMILES found for {zinc_id}"
            )
        
        # Create RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Invalid SMILES from ZINC: {smiles}"
            )
        
        # Standardize molecule
        mol = standardize_molecule(mol)
        
        # Add hydrogens if requested
        if add_hydrogens:
            mol = Chem.AddHs(mol)
        
        # Calculate properties
        properties = calculate_properties(mol)
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol)
        
        # Prepare response
        result = {
            "zinc_id": zinc_id,
            "smiles": canonical_smiles,
            "original_smiles": smiles,
            "properties": properties,
            "zinc_data": {
                "mwt": data.get("mwt"),
                "logp": data.get("logp"),
                "reactive": data.get("reactive"),
                "purchasability": data.get("purchasability"),
                "tranche_name": data.get("tranche_name")
            }
        }
        
        return types.TextContent(
            type="text",
            text=json.dumps(result, indent=2)
        )
        
    except requests.RequestException as e:
        return types.TextContent(
            type="text",
            text=f"Error fetching from ZINC: {str(e)}"
        )
    except Exception as e:
        logger.error(f"Error processing ZINC molecule {zinc_id}: {e}")
        return types.TextContent(
            type="text",
            text=f"Error processing molecule: {str(e)}"
        )


def fetch_zinc_batch(zinc_ids: str, add_hydrogens: bool = False) -> types.TextContent:
    """
    Fetch multiple molecules from ZINC database.
    
    Args:
        zinc_ids: Comma-separated list of ZINC IDs
        add_hydrogens: Whether to add explicit hydrogens
    
    Returns:
        TextContent with information for all molecules
    """
    try:
        ids = [id.strip() for id in zinc_ids.split(",")]
        results = []
        
        for zinc_id in ids:
            # Clean up ZINC ID
            zinc_id = zinc_id.upper()
            if not zinc_id.startswith("ZINC"):
                zinc_id = f"ZINC{zinc_id}"
            
            try:
                # Fetch molecule
                url = f"{ZINC_API_BASE}/substances/{zinc_id}"
                response = requests.get(url, timeout=10)
                
                if response.status_code == 404:
                    results.append({
                        "zinc_id": zinc_id,
                        "error": "Not found in ZINC database"
                    })
                    continue
                
                response.raise_for_status()
                data = response.json()
                
                smiles = data.get("smiles")
                if not smiles:
                    results.append({
                        "zinc_id": zinc_id,
                        "error": "No SMILES available"
                    })
                    continue
                
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    results.append({
                        "zinc_id": zinc_id,
                        "error": f"Invalid SMILES: {smiles}"
                    })
                    continue
                
                # Standardize
                mol = standardize_molecule(mol)
                
                # Add hydrogens if requested
                if add_hydrogens:
                    mol = Chem.AddHs(mol)
                
                # Calculate properties
                properties = calculate_properties(mol)
                
                results.append({
                    "zinc_id": zinc_id,
                    "smiles": Chem.MolToSmiles(mol),
                    "properties": properties,
                    "zinc_data": {
                        "mwt": data.get("mwt"),
                        "logp": data.get("logp"),
                        "purchasability": data.get("purchasability")
                    }
                })
                
            except Exception as e:
                results.append({
                    "zinc_id": zinc_id,
                    "error": str(e)
                })
        
        return types.TextContent(
            type="text",
            text=json.dumps({
                "count": len(results),
                "molecules": results
            }, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error in batch fetch: {e}")
        return types.TextContent(
            type="text",
            text=f"Error in batch fetch: {str(e)}"
        )


def search_zinc_similarity(smiles: str, similarity: float = 0.7, max_results: int = 10) -> types.TextContent:
    """
    Search ZINC database for similar molecules.
    
    Args:
        smiles: SMILES string of query molecule
        similarity: Tanimoto similarity threshold (0.0-1.0)
        max_results: Maximum number of results
    
    Returns:
        TextContent with similar molecules from ZINC
    """
    try:
        # Validate input molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Invalid SMILES string: {smiles}"
            )
        
        # Clean SMILES for query
        clean_smiles = Chem.MolToSmiles(mol)
        
        # ZINC similarity search endpoint
        url = f"{ZINC_API_BASE}/substances/search"
        params = {
            "q": clean_smiles,
            "type": "similarity",
            "threshold": similarity,
            "limit": max_results
        }
        
        response = requests.get(url, params=params, timeout=30)
        
        if response.status_code == 404:
            return types.TextContent(
                type="text",
                text="No similar molecules found in ZINC"
            )
        
        response.raise_for_status()
        data = response.json()
        
        results = []
        for item in data.get("results", []):
            zinc_id = item.get("zinc_id")
            mol_smiles = item.get("smiles")
            
            if mol_smiles:
                mol_obj = Chem.MolFromSmiles(mol_smiles)
                if mol_obj:
                    props = calculate_properties(mol_obj)
                    results.append({
                        "zinc_id": zinc_id,
                        "smiles": mol_smiles,
                        "similarity": item.get("similarity", similarity),
                        "properties": props
                    })
        
        return types.TextContent(
            type="text",
            text=json.dumps({
                "query_smiles": clean_smiles,
                "similarity_threshold": similarity,
                "count": len(results),
                "molecules": results
            }, indent=2)
        )
        
    except requests.RequestException as e:
        return types.TextContent(
            type="text",
            text=f"Error searching ZINC: {str(e)}"
        )
    except Exception as e:
        logger.error(f"Error in similarity search: {e}")
        return types.TextContent(
            type="text",
            text=f"Error in similarity search: {str(e)}"
        )


def search_zinc_substructure(smarts: str, max_results: int = 10) -> types.TextContent:
    """
    Search ZINC database for molecules containing a substructure.
    
    Args:
        smarts: SMARTS pattern for substructure search
        max_results: Maximum number of results
    
    Returns:
        TextContent with molecules containing the substructure
    """
    try:
        # Validate SMARTS pattern
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            return types.TextContent(
                type="text",
                text=f"Error: Invalid SMARTS pattern: {smarts}"
            )
        
        # ZINC substructure search endpoint
        url = f"{ZINC_API_BASE}/substances/search"
        params = {
            "q": smarts,
            "type": "substructure",
            "limit": max_results
        }
        
        response = requests.get(url, params=params, timeout=30)
        
        if response.status_code == 404:
            return types.TextContent(
                type="text",
                text="No molecules with this substructure found in ZINC"
            )
        
        response.raise_for_status()
        data = response.json()
        
        results = []
        for item in data.get("results", []):
            zinc_id = item.get("zinc_id")
            mol_smiles = item.get("smiles")
            
            if mol_smiles:
                mol_obj = Chem.MolFromSmiles(mol_smiles)
                if mol_obj:
                    # Verify substructure match
                    if mol_obj.HasSubstructMatch(pattern):
                        props = calculate_properties(mol_obj)
                        results.append({
                            "zinc_id": zinc_id,
                            "smiles": mol_smiles,
                            "properties": props
                        })
        
        return types.TextContent(
            type="text",
            text=json.dumps({
                "query_smarts": smarts,
                "count": len(results),
                "molecules": results
            }, indent=2)
        )
        
    except requests.RequestException as e:
        return types.TextContent(
            type="text",
            text=f"Error searching ZINC: {str(e)}"
        )
    except Exception as e:
        logger.error(f"Error in substructure search: {e}")
        return types.TextContent(
            type="text",
            text=f"Error in substructure search: {str(e)}"
        )


def standardize_molecule(mol: Chem.Mol) -> Chem.Mol:
    """
    Standardize a molecule using RDKit's standardization procedures.
    
    Args:
        mol: RDKit molecule object
    
    Returns:
        Standardized molecule
    """
    # Clean up the molecule
    mol = rdMolStandardize.Cleanup(mol)
    
    # Remove fragments (keep largest)
    mol = rdMolStandardize.FragmentParent(mol)
    
    # Neutralize charges where possible
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    
    # Sanitize
    Chem.SanitizeMol(mol)
    
    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    return mol


def calculate_properties(mol: Chem.Mol) -> Dict[str, Any]:
    """
    Calculate molecular properties for a molecule.
    
    Args:
        mol: RDKit molecule object
    
    Returns:
        Dictionary of molecular properties
    """
    return {
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_rings": Chem.GetSSSR(mol),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "num_chiral_centers": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    }


def prepare_zinc_molecule_for_editing(zinc_id: str, label_r_groups: bool = True) -> types.TextContent:
    """
    Fetch a ZINC molecule and prepare it for editing operations.
    
    Args:
        zinc_id: ZINC identifier
        label_r_groups: Whether to identify and label potential R-group sites
    
    Returns:
        TextContent with prepared molecule data including SMILES and edit sites
    """
    try:
        # Fetch the molecule
        result = fetch_zinc_molecule(zinc_id)
        
        # Parse the result
        data = json.loads(result.text)
        if "error" in data or "Error" in result.text:
            return result
        
        smiles = data["smiles"]
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            return types.TextContent(
                type="text",
                text=f"Error: Could not parse molecule {zinc_id}"
            )
        
        # Prepare for editing
        edit_info = {
            "zinc_id": zinc_id,
            "original_smiles": smiles,
            "properties": data["properties"]
        }
        
        if label_r_groups:
            # Identify potential edit sites
            # This is a simplified version - you might want more sophisticated detection
            
            # Find aromatic rings that could be substituted
            aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
            
            # Find terminal groups that could be replaced
            terminal_atoms = []
            for atom in mol.GetAtoms():
                if atom.GetDegree() == 1 and atom.GetAtomicNum() != 1:  # Not hydrogen
                    terminal_atoms.append(atom.GetIdx())
            
            edit_info["editable_sites"] = {
                "aromatic_positions": aromatic_atoms,
                "terminal_groups": terminal_atoms,
                "num_rings": Chem.GetSSSR(mol),
                "num_rotatable_bonds": data["properties"]["rotatable_bonds"]
            }
            
            # Generate a version with Bemis-Murcko scaffold
            from rdkit.Chem.Scaffolds import MurckoScaffold
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smiles = Chem.MolToSmiles(scaffold) if scaffold else None
            
            edit_info["scaffold"] = {
                "smiles": scaffold_smiles,
                "num_atoms": scaffold.GetNumAtoms() if scaffold else 0
            }
        
        edit_info["ready_for_editing"] = True
        
        return types.TextContent(
            type="text",
            text=json.dumps(edit_info, indent=2)
        )
        
    except Exception as e:
        logger.error(f"Error preparing ZINC molecule for editing: {e}")
        return types.TextContent(
            type="text",
            text=f"Error preparing molecule: {str(e)}"
        )
