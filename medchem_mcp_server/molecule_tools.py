"""
Molecule-related tools for the medchem MCP server.

This module contains all the molecular calculation and manipulation tools.
"""

import logging
import requests
import pandas

import mcp.types as types
from mcp.server.fastmcp import FastMCP

from medchem_mcp_server.utilities import pil_image_to_base64, img_byte_to_base64

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdRascalMCES

logger = logging.getLogger("medchem-mcp-server")


def get_molecule_descriptors_from_smiles(smiles: str) -> types.TextContent:
    """Get molecular descriptors from SMILES string"""
    
    m = Chem.MolFromSmiles(smiles)
    descriptors = [Descriptors.CalcMolDescriptors(m)]
    df = pandas.DataFrame(descriptors)
    print(df.head())
    table_as_string = df.to_string(index=False)
    return types.TextContent(
        type="text",
        text=table_as_string
    )


def search_substructure(smiles: str, substructure_smiles: str, use_chirality: bool = False) -> types.TextContent:
    """Search for substructure in a molecule"""
    m = Chem.MolFromSmiles(smiles)
    substructure = Chem.MolFromSmiles(substructure_smiles)
    match = m.HasSubstructMatch(substructure, useChirality=use_chirality)
    print("MATCH", match)
    is_match = None
    if match:
        is_match = True
        matches = m.GetSubstructMatches(substructure)
        return types.TextContent(
            type="text",
            text=f"Substructure found at: str {matches}",
        )
    else:
        is_match = False
        return types.TextContent(type="text", text="Substructure not found")


def get_smiles_from_name(chemical_name: str) -> types.TextContent:
    """Get SMILES string from name"""
    get_smile_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
    response = requests.get(get_smile_url)
    response_json = response.json()
    canonical_smiles = response_json["PropertyTable"]["Properties"][0].get("CanonicalSMILES", None)
    
    if canonical_smiles is None:
        canonical_smiles = response_json["PropertyTable"]["Properties"][0].get("ConnectivitySMILES", None)
    
    if canonical_smiles is None:
        return types.TextContent(type="text", text="Sorry No SMILES found")
    
    return types.TextContent(
        type="text",
        text=f"SMILES: {canonical_smiles}",
    )


def get_weight_from_smiles(smiles: str) -> types.TextContent:
    """Get molecular weight from SMILES string"""
    m = Chem.MolFromSmiles(smiles)
    weight = Descriptors.MolWt(m)
    return types.TextContent(
        type="text",
        text=f"Molecular Weight: {weight}",
    )


def find_maximum_common_substructure(smiles1: str, smiles2: str) -> types.TextContent:
    """Find maximum common substructure between two SMILES strings using RascalMCES"""
    m1 = Chem.MolFromSmiles(smiles1)
    m2 = Chem.MolFromSmiles(smiles2)
    mcs = rdRascalMCES.FindMCES(m1, m2)
    
    mcs_smiles = mcs.smartsString
    return types.TextContent(
        type="text",
        text=f"Maximum Common Substructure: {mcs_smiles}",
    )


def smiles_to_image(smiles: str) -> types.ImageContent:
    """Generate image from SMILES string"""
    m = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(m)
    print(dir(img))
    img_str = pil_image_to_base64(img)
    return types.ImageContent(
        type="image", data=img_str, mimeType="image/png"
    )
