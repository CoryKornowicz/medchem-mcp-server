"""
Molecule editing tools for medicinal chemistry.

This module provides safe, composable editing operations for molecules,
including R-group replacement, scaffold swapping, and substructure replacement.
Based on battle-tested patterns for LLM-guided molecule editing.
"""

import logging
import json
from typing import Optional, Dict, List, Any, Union, Tuple
import mcp.types as types
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import BRICS, Recap
from rdkit.Chem import Draw

logger = logging.getLogger("medchem-mcp-server")


# ============================================================================
# Core Loading and Preparation Functions
# ============================================================================

def load_and_prepare_molecule(smiles: str, sanitize: bool = True, 
                             remove_hs: bool = True) -> Optional[Chem.Mol]:
    """
    Load and prepare a molecule from SMILES with standard cleanup.
    
    Args:
        smiles: SMILES string
        sanitize: Whether to sanitize the molecule
        remove_hs: Whether to remove explicit hydrogens
    
    Returns:
        Prepared RDKit molecule or None if invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Optional cleanup and standardization
        mol = rdMolStandardize.Cleanup(mol)
        
        if sanitize:
            Chem.SanitizeMol(mol)
        
        if remove_hs:
            mol = Chem.RemoveHs(mol)
        
        # Assign stereochemistry
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        
        return mol
        
    except Exception as e:
        logger.error(f"Error preparing molecule: {e}")
        return None


# ============================================================================
# Scaffold and R-group Analysis
# ============================================================================

def get_bemis_murcko_scaffold(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """
    Extract Bemis-Murcko scaffold from a molecule.
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Scaffold molecule or None
    """
    try:
        return MurckoScaffold.GetScaffoldForMol(mol)
    except Exception as e:
        logger.error(f"Error extracting scaffold: {e}")
        return None


def label_r_groups(mol: Chem.Mol, core: Optional[Chem.Mol] = None) -> Optional[Chem.Mol]:
    """
    Replace core with dummy atoms labeled [*:1], [*:2], etc.
    This creates a representation of sidechains with attachment points.
    
    Args:
        mol: Original molecule
        core: Core/scaffold to replace (if None, uses Bemis-Murcko)
    
    Returns:
        Molecule with labeled R-groups or None
    """
    try:
        if core is None:
            core = get_bemis_murcko_scaffold(mol)
            if core is None:
                return None
        
        # Replace core with dummies
        sidechains = Chem.ReplaceCore(mol, core)
        if sidechains is None:
            return None
        
        # Add consistent labels [*:1], [*:2], ... to dummy atoms
        em = Chem.EditableMol(sidechains)
        sc = Chem.Mol(sidechains)
        
        label = 1
        for atom in sc.GetAtoms():
            if atom.GetAtomicNum() == 0:  # dummy atom
                atom.SetAtomMapNum(label)
                label += 1
        
        return sc
        
    except Exception as e:
        logger.error(f"Error labeling R-groups: {e}")
        return None


def create_templated_molecule(mol: Chem.Mol) -> Dict[str, Any]:
    """
    Create a templated version of a molecule with labeled attachment points.
    
    Args:
        mol: RDKit molecule
    
    Returns:
        Dictionary with template information
    """
    try:
        scaffold = get_bemis_murcko_scaffold(mol)
        if scaffold is None:
            return {"error": "Could not extract scaffold"}
        
        # Get sidechains with labels
        sidechains = label_r_groups(mol, scaffold)
        
        # Count R-groups
        r_groups = []
        if sidechains:
            for atom in sidechains.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0:
                    r_groups.append(atom.GetAtomMapNum())
        
        return {
            "original_smiles": Chem.MolToSmiles(mol),
            "scaffold_smiles": Chem.MolToSmiles(scaffold) if scaffold else None,
            "sidechain_smiles": Chem.MolToSmiles(sidechains) if sidechains else None,
            "num_r_groups": len(r_groups),
            "r_group_labels": sorted(r_groups)
        }
        
    except Exception as e:
        logger.error(f"Error creating template: {e}")
        return {"error": str(e)}


# ============================================================================
# Core Editing Operations
# ============================================================================

def replace_rgroup(mol: Chem.Mol, r: int, fragment: str) -> Optional[Chem.Mol]:
    """
    Replace dummy atom [*:r] in molecule with fragment containing [*:r].
    
    Args:
        mol: Molecule with labeled dummy atom
        r: R-group number
        fragment: SMILES of replacement fragment with [*:r]
    
    Returns:
        Edited molecule or None
    """
    try:
        # Find attachment site
        query = Chem.MolFromSmarts(f"[*:{r}]")
        if query is None:
            return None
        
        # Parse replacement fragment
        replacement = Chem.MolFromSmiles(fragment)
        if replacement is None:
            return None
        
        # Perform replacement
        results = AllChem.ReplaceSubstructs(mol, query, replacement, replaceAll=False)
        
        if results and len(results) > 0:
            result_mol = results[0]
            Chem.SanitizeMol(result_mol)
            return result_mol
        
        return None
        
    except Exception as e:
        logger.error(f"Error in R-group replacement: {e}")
        return None


def substructure_replace(mol: Chem.Mol, query_smarts: str, 
                        replacement_smiles: str, 
                        replace_all: bool = False) -> List[Chem.Mol]:
    """
    Replace substructure matches with new fragment, preserving attachments.
    
    Args:
        mol: Target molecule
        query_smarts: SMARTS pattern to find (can include [*:n] labels)
        replacement_smiles: Replacement SMILES (should match [*:n] labels)
        replace_all: Replace all occurrences or just first
    
    Returns:
        List of resulting molecules
    """
    try:
        query = Chem.MolFromSmarts(query_smarts)
        if query is None:
            return []
        
        replacement = Chem.MolFromSmiles(replacement_smiles)
        if replacement is None:
            return []
        
        results = AllChem.ReplaceSubstructs(mol, query, replacement, replaceAll=replace_all)
        
        # Sanitize results
        valid_results = []
        for res_mol in results:
            try:
                Chem.SanitizeMol(res_mol)
                valid_results.append(res_mol)
            except:
                continue
        
        return valid_results
        
    except Exception as e:
        logger.error(f"Error in substructure replacement: {e}")
        return []


def scaffold_swap(mol: Chem.Mol, old_core_smarts: str, 
                 new_core_smiles: str) -> Optional[Chem.Mol]:
    """
    Swap scaffold/core while preserving substituents via labeled atoms.
    
    Args:
        mol: Target molecule
        old_core_smarts: SMARTS of current core with [*:n] labels
        new_core_smiles: SMILES of new core with matching [*:n] labels
    
    Returns:
        Modified molecule or None
    """
    try:
        query = Chem.MolFromSmarts(old_core_smarts)
        if query is None:
            return None
        
        new_core = Chem.MolFromSmiles(new_core_smiles)
        if new_core is None:
            return None
        
        results = AllChem.ReplaceSubstructs(mol, query, new_core, replaceAll=False)
        
        if results and len(results) > 0:
            result_mol = results[0]
            Chem.SanitizeMol(result_mol)
            return result_mol
        
        return None
        
    except Exception as e:
        logger.error(f"Error in scaffold swap: {e}")
        return None


# ============================================================================
# Edit Executor (JSON-based edit plans)
# ============================================================================

def apply_edit_plan(mol: Chem.Mol, edits_json: str) -> Tuple[Optional[Chem.Mol], List[str]]:
    """
    Apply a JSON-encoded list of edit operations to a molecule.
    
    Args:
        mol: Starting molecule
        edits_json: JSON string with list of edit operations
    
    Returns:
        Tuple of (edited molecule, list of status messages)
    """
    messages = []
    
    try:
        edits = json.loads(edits_json)
        if not isinstance(edits, list):
            return None, ["Error: edits must be a JSON array"]
        
        current_mol = Chem.Mol(mol)  # Make a copy
        
        for i, op in enumerate(edits):
            op_type = op.get("op")
            
            if op_type == "rgroup_replace":
                r = int(op.get("r", 0))
                fragment = op.get("fragment", "")
                
                result = replace_rgroup(current_mol, r, fragment)
                if result:
                    current_mol = result
                    messages.append(f"Step {i+1}: Replaced R{r} with {fragment}")
                else:
                    messages.append(f"Step {i+1}: Failed to replace R{r}")
                    
            elif op_type == "substructure_replace":
                query = op.get("query_smarts", "")
                replacement = op.get("replacement_smiles", "")
                replace_all = op.get("replace_all", False)
                
                results = substructure_replace(current_mol, query, replacement, replace_all)
                if results:
                    current_mol = results[0]
                    messages.append(f"Step {i+1}: Replaced substructure {query}")
                else:
                    messages.append(f"Step {i+1}: Failed to replace substructure")
                    
            elif op_type == "scaffold_swap":
                old_core = op.get("old_core_smarts", "")
                new_core = op.get("new_core_smiles", "")
                
                result = scaffold_swap(current_mol, old_core, new_core)
                if result:
                    current_mol = result
                    messages.append(f"Step {i+1}: Swapped scaffold")
                else:
                    messages.append(f"Step {i+1}: Failed to swap scaffold")
                    
            else:
                messages.append(f"Step {i+1}: Unknown operation '{op_type}'")
            
            # Validate after each step
            try:
                Chem.SanitizeMol(current_mol)
                Chem.AssignStereochemistry(current_mol, cleanIt=True, force=True)
            except Exception as e:
                messages.append(f"Step {i+1}: Sanitization warning: {e}")
        
        return current_mol, messages
        
    except json.JSONDecodeError:
        return None, ["Error: Invalid JSON format"]
    except Exception as e:
        return None, [f"Error applying edits: {str(e)}"]


# ============================================================================
# Validation and Filtering
# ============================================================================

def validate_molecule(mol: Chem.Mol, filters: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Validate molecule and apply property filters.
    
    Args:
        mol: Molecule to validate
        filters: Dictionary of property filters (e.g., {"mw_max": 500})
    
    Returns:
        Validation results dictionary
    """
    if filters is None:
        filters = {
            "mw_max": 700,
            "mw_min": 150,
            "logp_max": 6,
            "logp_min": -2,
            "hbd_max": 5,
            "hba_max": 10,
            "rotatable_bonds_max": 10,
            "tpsa_max": 140
        }
    
    results = {"valid": True, "warnings": [], "properties": {}}
    
    try:
        # Calculate properties
        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "tpsa": Descriptors.TPSA(mol),
            "qed": Descriptors.qed(mol),
            "sa_score": None  # Would need external SA score calculator
        }
        results["properties"] = props
        
        # Check filters
        if "mw_max" in filters and props["mw"] > filters["mw_max"]:
            results["warnings"].append(f"MW {props['mw']:.1f} exceeds limit {filters['mw_max']}")
            results["valid"] = False
            
        if "mw_min" in filters and props["mw"] < filters["mw_min"]:
            results["warnings"].append(f"MW {props['mw']:.1f} below minimum {filters['mw_min']}")
            results["valid"] = False
            
        if "logp_max" in filters and props["logp"] > filters["logp_max"]:
            results["warnings"].append(f"LogP {props['logp']:.2f} exceeds limit {filters['logp_max']}")
            results["valid"] = False
            
        if "hbd_max" in filters and props["hbd"] > filters["hbd_max"]:
            results["warnings"].append(f"HBD {props['hbd']} exceeds limit {filters['hbd_max']}")
            results["valid"] = False
            
        if "hba_max" in filters and props["hba"] > filters["hba_max"]:
            results["warnings"].append(f"HBA {props['hba']} exceeds limit {filters['hba_max']}")
            results["valid"] = False
            
        if "rotatable_bonds_max" in filters and props["rotatable_bonds"] > filters["rotatable_bonds_max"]:
            results["warnings"].append(f"Rotatable bonds {props['rotatable_bonds']} exceeds limit")
            results["valid"] = False
            
        if "tpsa_max" in filters and props["tpsa"] > filters["tpsa_max"]:
            results["warnings"].append(f"TPSA {props['tpsa']:.1f} exceeds limit {filters['tpsa_max']}")
            results["valid"] = False
        
        # Lipinski's Rule of 5
        ro5_violations = 0
        if props["mw"] > 500:
            ro5_violations += 1
        if props["logp"] > 5:
            ro5_violations += 1
        if props["hbd"] > 5:
            ro5_violations += 1
        if props["hba"] > 10:
            ro5_violations += 1
            
        if ro5_violations > 1:
            results["warnings"].append(f"Lipinski Ro5: {ro5_violations} violations")
        
        results["ro5_violations"] = ro5_violations
        
    except Exception as e:
        results["valid"] = False
        results["warnings"].append(f"Validation error: {str(e)}")
    
    return results


# ============================================================================
# Advanced Editing Tools
# ============================================================================

def fragment_molecule(mol: Chem.Mol, method: str = "brics") -> List[str]:
    """
    Fragment molecule using BRICS or RECAP.
    
    Args:
        mol: Molecule to fragment
        method: "brics" or "recap"
    
    Returns:
        List of fragment SMILES
    """
    try:
        if method.lower() == "brics":
            fragments = BRICS.BRICSDecompose(mol)
        elif method.lower() == "recap":
            hierarchy = Recap.RecapDecompose(mol)
            fragments = [Chem.MolToSmiles(n.mol) for n in hierarchy.GetLeaves()]
        else:
            return []
        
        return list(fragments)
        
    except Exception as e:
        logger.error(f"Error fragmenting molecule: {e}")
        return []


def enumerate_r_groups(core_smiles: str, r_group_library: List[str], 
                       r_position: int = 1) -> List[str]:
    """
    Enumerate molecules by attaching R-groups from a library.
    
    Args:
        core_smiles: Core with [*:n] attachment point
        r_group_library: List of fragment SMILES with [*:n]
        r_position: Which R-group position to enumerate
    
    Returns:
        List of enumerated SMILES
    """
    results = []
    
    try:
        core = Chem.MolFromSmiles(core_smiles)
        if core is None:
            return []
        
        for fragment in r_group_library:
            mol_copy = Chem.Mol(core)
            result = replace_rgroup(mol_copy, r_position, fragment)
            
            if result:
                smiles = Chem.MolToSmiles(result)
                results.append(smiles)
                
    except Exception as e:
        logger.error(f"Error in R-group enumeration: {e}")
    
    return results


def apply_transformation(mol: Chem.Mol, smirks: str) -> Optional[Chem.Mol]:
    """
    Apply a SMIRKS transformation to a molecule.
    
    Args:
        mol: Target molecule
        smirks: SMIRKS reaction pattern
    
    Returns:
        Transformed molecule or None
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smirks)
        if rxn is None:
            return None
        
        products = rxn.RunReactants((mol,))
        
        if products and len(products) > 0:
            product = products[0][0]
            Chem.SanitizeMol(product)
            return product
        
        return None
        
    except Exception as e:
        logger.error(f"Error applying transformation: {e}")
        return None


# ============================================================================
# Common Fragment Libraries
# ============================================================================

# Common R-groups for medicinal chemistry
COMMON_R_GROUPS = {
    "alkyl": [
        "[*:1]C", "[*:1]CC", "[*:1]CCC", "[*:1]C(C)C", "[*:1]C(C)(C)C",
        "[*:1]C1CC1", "[*:1]C1CCC1", "[*:1]C1CCCC1"
    ],
    "aromatic": [
        "[*:1]c1ccccc1", "[*:1]c1ccc(F)cc1", "[*:1]c1ccc(Cl)cc1",
        "[*:1]c1ccc(C)cc1", "[*:1]c1ccc(OC)cc1", "[*:1]c1ccc(CF3)cc1",
        "[*:1]c1ccncc1", "[*:1]c1ncccn1", "[*:1]c1cccnc1"
    ],
    "heteroaromatic": [
        "[*:1]c1nccs1", "[*:1]c1ncco1", "[*:1]c1[nH]ccc1",
        "[*:1]c1n[nH]cc1", "[*:1]c1noc(C)c1"
    ],
    "polar": [
        "[*:1]O", "[*:1]OC", "[*:1]OCC", "[*:1]N", "[*:1]NC",
        "[*:1]N(C)C", "[*:1]S", "[*:1]SC", "[*:1]S(=O)(=O)C"
    ],
    "halogen": [
        "[*:1]F", "[*:1]Cl", "[*:1]Br", "[*:1]I"
    ]
}

# Common scaffolds for med chem
COMMON_SCAFFOLDS = {
    "phenyl": "[*:1]c1ccc([*:2])cc1",
    "pyridyl": "[*:1]c1ncc([*:2])cc1",
    "pyrimidyl": "[*:1]c1ncc([*:2])cn1",
    "thiophene": "[*:1]c1sc([*:2])cc1",
    "furan": "[*:1]c1oc([*:2])cc1",
    "pyrrole": "[*:1]c1[nH]c([*:2])cc1",
    "imidazole": "[*:1]c1nc([*:2])[nH]c1",
    "oxazole": "[*:1]c1nc([*:2])oc1",
    "thiazole": "[*:1]c1nc([*:2])sc1",
    "indole": "[*:1]c1c[nH]c2c([*:2])cccc12",
    "benzimidazole": "[*:1]c1nc2c([*:2])cccc2[nH]1",
    "quinoline": "[*:1]c1nc2c([*:2])cccc2c1"
}

# Common transformations (SMIRKS)
COMMON_TRANSFORMATIONS = {
    "fluorinate_para": "[c:1][cH:2][c:3]>>[c:1][c:2](F)[c:3]",
    "methylate_nitrogen": "[NH:1]>>[N:1](C)",
    "acetylate_amine": "[NH2:1]>>[NH:1](C(=O)C)",
    "oxidize_sulfur": "[#16:1]>>[#16:1](=O)",
    "reduce_nitro": "[N+:1](=O)[O-]>>[NH2:1]",
    "brominate_aromatic": "[cH:1]>>[c:1](Br)",
    "form_amide": "[C:1](=O)O.[N:2]>>[C:1](=O)[N:2]"
}
