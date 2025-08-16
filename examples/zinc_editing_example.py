#!/usr/bin/env python3
"""
Example script demonstrating ZINC molecule fetching and editing capabilities.

This example shows how to:
1. Fetch molecules from ZINC database
2. Apply various editing operations
3. Create analog series
4. Perform scaffold hopping
"""

import json
from medchem_mcp_server.zinc_tooling import (
    fetch_zinc_molecule,
    search_zinc_similarity
)
from medchem_mcp_server.zinc_editing_tools import (
    fetch_and_edit_zinc_molecule,
    create_zinc_analog_series,
    scaffold_hop_zinc_molecule,
    apply_medchem_transforms
)
from medchem_mcp_server.molecule_editing import (
    load_and_prepare_molecule,
    validate_molecule
)


def example_basic_fetch():
    """Example: Fetch a molecule from ZINC"""
    print("\n=== Example 1: Fetching a ZINC Molecule ===")
    
    # Fetch a known drug-like molecule
    result = fetch_zinc_molecule("ZINC000000000001")
    data = json.loads(result.text)
    
    print(f"ZINC ID: {data['zinc_id']}")
    print(f"SMILES: {data['smiles']}")
    print(f"MW: {data['properties']['molecular_weight']}")
    print(f"LogP: {data['properties']['logp']}")


def example_edit_zinc_molecule():
    """Example: Fetch and edit a ZINC molecule"""
    print("\n=== Example 2: Editing a ZINC Molecule ===")
    
    # Define edit operations
    edits = [
        {
            "op": "substructure_replace",
            "query_smarts": "c1ccccc1",  # Replace benzene
            "replacement_smiles": "c1ncccn1",  # With pyrimidine
            "replace_all": False
        }
    ]
    
    # Apply edits to a ZINC molecule
    result = fetch_and_edit_zinc_molecule(
        zinc_id="ZINC000000000001",
        edits_json=json.dumps(edits),
        session_id="example_session"
    )
    
    data = json.loads(result.text)
    print(f"Original SMILES: {data['original_smiles']}")
    print(f"Edited SMILES: {data['edited_smiles']}")
    print(f"Validation: {data['validation']['valid']}")
    if data['validation'].get('warnings'):
        print(f"Warnings: {data['validation']['warnings']}")


def example_analog_series():
    """Example: Create an analog series"""
    print("\n=== Example 3: Creating Analog Series ===")
    
    # Create aromatic analogs
    result = create_zinc_analog_series(
        zinc_id="ZINC000000000001",
        r_group_type="aromatic",
        r_position=1,
        session_id="example_session"
    )
    
    data = json.loads(result.text)
    print(f"Parent molecule: {data['parent_zinc_id']}")
    print(f"Number of analogs: {data['num_analogs']}")
    
    # Show first 3 valid analogs
    valid_analogs = [a for a in data['analogs'] if a['valid']][:3]
    for i, analog in enumerate(valid_analogs, 1):
        print(f"\nAnalog {i}:")
        print(f"  SMILES: {analog['smiles']}")
        print(f"  MW: {analog['properties']['molecular_weight']:.1f}")
        print(f"  LogP: {analog['properties']['logp']:.2f}")


def example_scaffold_hopping():
    """Example: Scaffold hopping"""
    print("\n=== Example 4: Scaffold Hopping ===")
    
    # Replace phenyl scaffold with pyridyl
    result = scaffold_hop_zinc_molecule(
        zinc_id="ZINC000000000001",
        new_scaffold_type="pyridyl",
        session_id="example_session"
    )
    
    data = json.loads(result.text)
    
    if "error" not in data and "Error" not in str(data):
        print(f"Original scaffold: {data.get('original_scaffold', 'N/A')}")
        print(f"New scaffold type: {data['new_scaffold_type']}")
        print(f"Hopped molecule: {data.get('hopped_smiles', 'N/A')}")
        if data.get('validation'):
            print(f"Valid: {data['validation']['valid']}")


def example_medchem_transforms():
    """Example: Apply medicinal chemistry transformations"""
    print("\n=== Example 5: MedChem Transformations ===")
    
    # Apply common transformations
    result = apply_medchem_transforms(
        zinc_id="ZINC000000000001",
        transform_names="fluorinate_para,methylate_nitrogen",
        session_id="example_session"
    )
    
    data = json.loads(result.text)
    print(f"Original molecule: {data['original_smiles']}")
    
    for transform in data['transformations']:
        if 'error' not in transform:
            print(f"\nTransformation: {transform['transform']}")
            print(f"  SMIRKS: {transform['smirks']}")
            print(f"  Result: {transform.get('smiles', 'N/A')}")
            print(f"  Valid: {transform.get('valid', False)}")


def example_similarity_search():
    """Example: Search for similar molecules in ZINC"""
    print("\n=== Example 6: Similarity Search ===")
    
    # Search for molecules similar to aspirin
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    
    result = search_zinc_similarity(
        smiles=aspirin_smiles,
        similarity=0.7,
        max_results=5
    )
    
    data = json.loads(result.text)
    print(f"Query molecule: {data['query_smiles']}")
    print(f"Similarity threshold: {data['similarity_threshold']}")
    print(f"Found {data['count']} similar molecules")
    
    for mol in data.get('molecules', [])[:3]:
        print(f"\n  ZINC ID: {mol['zinc_id']}")
        print(f"  Similarity: {mol['similarity']:.2f}")
        print(f"  MW: {mol['properties']['molecular_weight']:.1f}")


def example_complex_edit_plan():
    """Example: Complex multi-step editing"""
    print("\n=== Example 7: Complex Edit Plan ===")
    
    # Multi-step edit plan
    edits = [
        {
            "op": "rgroup_replace",
            "r": 1,
            "fragment": "c1ccc(F)cc1[*:1]"  # Add fluorobenzene
        },
        {
            "op": "substructure_replace",
            "query_smarts": "[NH2:1]",  # Find primary amine
            "replacement_smiles": "[NH:1](C)",  # Replace with methylamine
            "replace_all": True
        },
        {
            "op": "scaffold_swap",
            "old_core_smarts": "[*:1]c1ccccc1[*:2]",  # Phenyl scaffold
            "new_core_smiles": "[*:1]c1ncccn1[*:2]"  # Pyrimidine scaffold
        }
    ]
    
    result = fetch_and_edit_zinc_molecule(
        zinc_id="ZINC000000000001",
        edits_json=json.dumps(edits),
        session_id="example_session"
    )
    
    data = json.loads(result.text)
    print(f"Edit messages:")
    for msg in data.get('edit_messages', []):
        print(f"  - {msg}")
    
    if 'edited_smiles' in data:
        print(f"\nFinal molecule: {data['edited_smiles']}")
        
        # Validate the result
        mol = load_and_prepare_molecule(data['edited_smiles'])
        if mol:
            validation = validate_molecule(mol)
            print(f"Lipinski Ro5 violations: {validation.get('ro5_violations', 0)}")


if __name__ == "__main__":
    print("ZINC Molecule Fetching and Editing Examples")
    print("=" * 50)
    
    # Note: These examples assume ZINC API is accessible
    # Some may fail if ZINC IDs don't exist or API is unavailable
    
    try:
        example_basic_fetch()
    except Exception as e:
        print(f"Error in basic fetch: {e}")
    
    try:
        example_edit_zinc_molecule()
    except Exception as e:
        print(f"Error in editing: {e}")
    
    try:
        example_analog_series()
    except Exception as e:
        print(f"Error in analog series: {e}")
    
    try:
        example_scaffold_hopping()
    except Exception as e:
        print(f"Error in scaffold hopping: {e}")
    
    try:
        example_medchem_transforms()
    except Exception as e:
        print(f"Error in transformations: {e}")
    
    try:
        example_similarity_search()
    except Exception as e:
        print(f"Error in similarity search: {e}")
    
    try:
        example_complex_edit_plan()
    except Exception as e:
        print(f"Error in complex editing: {e}")
    
    print("\n" + "=" * 50)
    print("Examples completed!")
