#!/usr/bin/env python3
"""
Example workflow: 1LEP PDB structure docking with AutoDock Vina

This script demonstrates the complete workflow for molecular docking using the 
medchem MCP server tools:

1. Fetch 1LEP PDB structure (Human carbonic anhydrase II with acetazolamide)
2. Separate protein and ligand components
3. Prepare protein for docking (fix structure, add hydrogens, convert to PDBQT)
4. Prepare ligand for docking (convert to SDF, then PDBQT)
5. Perform molecular docking with AutoDock Vina
6. Score the binding pose

1LEP contains:
- Protein: Human carbonic anhydrase II (HCAII)
- Ligand: Acetazolamide (AZM) - a carbonic anhydrase inhibitor
"""

import sys
import os

# Add the parent directory to the path so we can import the server modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from medchem_mcp_server.protein_tools import fetch_pdb_file
from medchem_mcp_server.session_state import session_manager
from medchem_mcp_server.session_tools import (
    prepare_protein_for_docking,
    prepare_ligand_for_docking,
    dock_ligand_with_vina,
    score_ligand_with_vina,
    list_stored_files
)
from medchem_mcp_server.utilities import separate_protein_ligand_pdb


def print_section(title: str):
    """Print a formatted section header"""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def print_step(step: str):
    """Print a formatted step"""
    print(f"\n🔹 {step}")


def main():
    """Main workflow for 1LEP docking example"""
    
    print("🧬 1LEP Molecular Docking Workflow with AutoDock Vina")
    print("Structure: Human Carbonic Anhydrase II + Acetazolamide")
    
    try:
        # Step 1: Create a session
        print_section("Step 1: Session Setup")
        print_step("Creating new session...")
        session = session_manager.create_session()
        session_id = session.session_id
        print(f"✅ Session created: {session_id}")
        
        # Step 2: Fetch 1LEP PDB structure
        print_section("Step 2: Fetch PDB Structure")
        print_step("Downloading 1LEP from RCSB PDB...")
        fetch_result = fetch_pdb_file(
            pdb_id="1lep",
            session_id=session_id,
            save_to_session=True,
            tags="carbonic_anhydrase,acetazolamide,example,docking"
        )
        
        if "Error" in fetch_result.text:
            print(f"❌ Failed to fetch PDB: {fetch_result.text}")
            return False
        
        print("✅ 1LEP PDB structure downloaded successfully")
        
        # Step 3: List files to find the PDB file ID
        print_section("Step 3: Locate PDB File")
        print_step("Finding the downloaded PDB file...")
        files_result = list_stored_files(session_id)
        
        # Parse the file list to find the 1LEP PDB file
        pdb_file_id = None
        for line in files_result.text.split('\n'):
            if '1lep.pdb' in line and 'File ID:' in line:
                pdb_file_id = line.split('File ID: ')[1].split()[0]
                break
        
        if not pdb_file_id:
            print("❌ Could not find the PDB file ID")
            return False
            
        print(f"✅ Found PDB file ID: {pdb_file_id}")
        
        # Step 4: Separate protein and ligand components
        print_section("Step 4: Separate Protein and Ligand")
        print_step("Extracting protein and ligand components...")
        
        # Get the PDB file path
        pdb_file = session_manager.get_file(session_id, pdb_file_id)
        if not pdb_file:
            print("❌ Could not retrieve PDB file")
            return False
        
        # Separate protein and ligand
        protein_path, ligand_paths = separate_protein_ligand_pdb(
            pdb_file.file_path, 
            output_dir=session.tmp_dir
        )
        
        if not protein_path or not ligand_paths:
            print("❌ Failed to separate protein and ligand")
            return False
        
        print(f"✅ Protein extracted: {os.path.basename(protein_path)}")
        print(f"✅ Ligand(s) extracted: {[os.path.basename(p) for p in ligand_paths]}")
        
        # Save the separated files to the session
        # Save protein
        with open(protein_path, 'r') as f:
            protein_content = f.read()
        
        protein_file_id = session_manager.save_file(
            session_id=session_id,
            content=protein_content,
            filename="1lep_protein.pdb",
            file_type="pdb",
            tags=["protein", "1lep", "carbonic_anhydrase", "separated"]
        )
        
        # Save ligand (use the first one if multiple)
        ligand_path = ligand_paths[0]
        with open(ligand_path, 'r') as f:
            ligand_content = f.read()
        
        ligand_file_id = session_manager.save_file(
            session_id=session_id,
            content=ligand_content,
            filename="1lep_ligand.pdb",
            file_type="pdb",
            tags=["ligand", "1lep", "acetazolamide", "separated"]
        )
        
        print(f"✅ Protein file ID: {protein_file_id}")
        print(f"✅ Ligand file ID: {ligand_file_id}")
        
        # Step 5: Prepare protein for docking
        print_section("Step 5: Prepare Protein for Docking")
        print_step("Processing protein: fixing structure, adding hydrogens, converting to PDBQT...")
        
        protein_prep_result = prepare_protein_for_docking(
            session_id=session_id,
            file_id=protein_file_id,
            fix_structure=True,
            add_hydrogens=True,
            ph=7.0
        )
        
        print(protein_prep_result.text)
        
        # Extract protein PDBQT file ID from the result
        protein_pdbqt_id = None
        for line in protein_prep_result.text.split('\n'):
            if 'PDBQT file ID:' in line:
                protein_pdbqt_id = line.split('PDBQT file ID: ')[1].strip()
                break
        
        if not protein_pdbqt_id:
            print("❌ Failed to prepare protein PDBQT")
            return False
        
        print(f"✅ Protein PDBQT file ID: {protein_pdbqt_id}")
        
        # Step 6: Prepare ligand for docking
        print_section("Step 6: Prepare Ligand for Docking")
        print_step("Converting ligand to SDF, then to PDBQT...")
        
        ligand_prep_result = prepare_ligand_for_docking(
            session_id=session_id,
            file_id=ligand_file_id,
            name="acetazolamide"
        )
        
        print(ligand_prep_result.text)
        
        # Extract ligand PDBQT file ID from the result
        ligand_pdbqt_id = None
        for line in ligand_prep_result.text.split('\n'):
            if 'PDBQT file ID:' in line:
                ligand_pdbqt_id = line.split('PDBQT file ID: ')[1].strip()
                break
        
        if not ligand_pdbqt_id:
            print("❌ Failed to prepare ligand PDBQT")
            return False
        
        print(f"✅ Ligand PDBQT file ID: {ligand_pdbqt_id}")
        
        # Step 7: Perform molecular docking
        print_section("Step 7: Molecular Docking with AutoDock Vina")
        print_step("Running molecular docking (this may take a few minutes)...")
        
        docking_result = dock_ligand_with_vina(
            session_id=session_id,
            protein_file_id=protein_pdbqt_id,
            ligand_file_id=ligand_pdbqt_id,
            center=None,  # Will be calculated from ligand center of mass
            box_size=[20.0, 20.0, 20.0],  # 20 Å box
            exhaustiveness=8,
            num_modes=9,
            energy_range=3.0
        )
        
        print(docking_result.text)
        
        # Step 8: Score the original ligand pose
        print_section("Step 8: Score Original Ligand Pose")
        print_step("Scoring the crystal structure ligand pose...")
        
        scoring_result = score_ligand_with_vina(
            session_id=session_id,
            protein_file_id=protein_pdbqt_id,
            ligand_file_id=ligand_pdbqt_id,
            center=None,  # Will be calculated from ligand center of mass
            box_size=[20.0, 20.0, 20.0]
        )
        
        print(scoring_result.text)
        
        # Step 9: Summary
        print_section("Step 9: Workflow Summary")
        print("🎯 1LEP Docking Workflow Completed Successfully!")
        print("\nFiles created in this session:")
        
        final_files = list_stored_files(session_id)
        print(final_files.text)
        
        print(f"\n💡 Session ID for future reference: {session_id}")
        print("\nYou can now:")
        print("  • Analyze the docked poses")
        print("  • Compare binding affinities")
        print("  • Visualize the results")
        print("  • Run additional docking experiments")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error in workflow: {str(e)}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Starting 1LEP docking workflow...")
    success = main()
    
    if success:
        print("\n🎉 Workflow completed successfully!")
        sys.exit(0)
    else:
        print("\n💥 Workflow failed!")
        sys.exit(1)
