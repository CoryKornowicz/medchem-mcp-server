#!/usr/bin/env python3
"""
Test script for session management functionality
"""

from medchem_mcp_server.session_state import session_manager
from medchem_mcp_server.session_tools import (
    create_session,
    store_molecule,
    store_protein,
    list_stored_molecules,
    list_stored_proteins,
    get_session_summary
)

def test_session_management():
    """Test basic session management operations"""
    
    print("Testing Session Management System")
    print("=" * 50)
    
    # Create a new session
    result = create_session()
    print("\n1. Creating new session:")
    print(result.text)
    
    # Extract session ID from the result
    session_id = result.text.split("Session ID: ")[1].split("\n")[0]
    print(f"\nExtracted Session ID: {session_id}")
    
    # Store a molecule
    print("\n2. Storing a molecule (aspirin):")
    result = store_molecule(
        smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        name="Aspirin",
        session_id=session_id,
        tags="nsaid,pain-relief"
    )
    print(result.text)
    
    # Store another molecule
    print("\n3. Storing another molecule (caffeine):")
    result = store_molecule(
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        name="Caffeine",
        session_id=session_id,
        tags="stimulant"
    )
    print(result.text)
    
    # Store a protein (without fetching, just metadata)
    print("\n4. Storing protein metadata:")
    result = store_protein(
        pdb_id="1hsg",
        name="HIV-1 Protease",
        session_id=session_id,
        tags="protease,hiv"
    )
    print(result.text)
    
    # List stored molecules
    print("\n5. Listing stored molecules:")
    result = list_stored_molecules(session_id=session_id)
    print(result.text)
    
    # List stored proteins
    print("\n6. Listing stored proteins:")
    result = list_stored_proteins(session_id=session_id)
    print(result.text)
    
    # Get session summary
    print("\n7. Session summary:")
    result = get_session_summary(session_id=session_id)
    print(result.text)
    
    # Test file storage via fetch_pdb_file
    print("\n8. Testing fetch_pdb_file with session storage:")
    from medchem_mcp_server.protein_tools import fetch_pdb_file
    result = fetch_pdb_file("1crn", session_id=session_id, tags="small,crambin")
    # Only print first 500 chars and session info
    content = result.text
    if len(content) > 500:
        # Find the session info part
        if "[Saved to session" in content:
            session_info_start = content.index("[Saved to session")
            print(content[:500] + "\n...\n" + content[session_info_start:])
        else:
            print(content[:500] + "\n...")
    else:
        print(content)
    
    # Final session summary
    print("\n9. Final session summary:")
    result = get_session_summary(session_id=session_id)
    print(result.text)
    
    print("\n" + "=" * 50)
    print("Test completed successfully!")

if __name__ == "__main__":
    test_session_management()
