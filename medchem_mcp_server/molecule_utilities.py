"""
Molecule utilities for converting between formats and extracting chemical information.
"""

import os
import tempfile
import logging
from typing import Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from io import StringIO
import pypdb
from prody import parsePDB, writePDBStream
from openmm.app import PDBFile


logger = logging.getLogger("medchem-mcp-server")


def pdb_to_sdf(pdb_content: str, res_name: Optional[str] = None, 
               remove_waters: bool = True, remove_ions: bool = True) -> Optional[str]:
    """
    Convert PDB content to SDF format using ProDy and RDKit with proper bond order assignment.
    
    Args:
        pdb_content: PDB file content as string
        res_name: Residue name for bond order assignment (if known)
        remove_waters: Remove water molecules (HOH, WAT)
        remove_ions: Remove common ions
        
    Returns:
        SDF content as string, or None if conversion failed
    """
    try:
        # If we have a residue name, use ProDy and the helper functions for proper bond orders
        if res_name:
            try:
                # Parse with ProDy
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
                    tmp_pdb.write(pdb_content)
                    tmp_pdb_path = tmp_pdb.name
                
                try:
                    ligand = parsePDB(tmp_pdb_path)
                    if ligand:
                        # Process with proper bond order assignment
                        mols = process_ligand_stereochemistry(ligand, res_name)
                        if mols:
                            # Use the first conformer
                            mol = mols[0]
                            # Convert to SDF
                            sdf_content = Chem.MolToMolBlock(mol)
                            return sdf_content
                finally:
                    if os.path.exists(tmp_pdb_path):
                        os.unlink(tmp_pdb_path)
            except Exception as e:
                logger.warning(f"Could not process with ProDy/pypdb, falling back to basic conversion: {e}")
        
        # Fallback to basic RDKit conversion
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
            tmp_pdb.write(pdb_content)
            tmp_pdb_path = tmp_pdb.name
        
        try:
            # Read PDB file with RDKit
            mol = Chem.MolFromPDBFile(tmp_pdb_path, removeHs=False, sanitize=False)
            
            if mol is None:
                logger.error("Failed to parse PDB file with RDKit")
                return None
            
            # Remove waters and ions if requested
            if remove_waters or remove_ions:
                editable_mol = Chem.EditableMol(mol)
                atoms_to_remove = []
                for atom in mol.GetAtoms():
                    residue_info = atom.GetPDBResidueInfo()
                    if residue_info:
                        res_name_atom = residue_info.GetResidueName().strip()
                        if remove_waters and res_name_atom in ['HOH', 'WAT', 'H2O']:
                            atoms_to_remove.append(atom.GetIdx())
                        elif remove_ions and res_name_atom in ['NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE']:
                            atoms_to_remove.append(atom.GetIdx())
                
                # Remove atoms in reverse order
                for idx in sorted(atoms_to_remove, reverse=True):
                    editable_mol.RemoveAtom(idx)
                mol = editable_mol.GetMol()
            
            # Try to sanitize
            try:
                Chem.SanitizeMol(mol)
            except:
                logger.warning("Could not fully sanitize molecule")
            
            # Convert to SDF
            sdf_content = Chem.MolToMolBlock(mol)
            return sdf_content
            
        finally:
            # Clean up temp file
            if os.path.exists(tmp_pdb_path):
                os.unlink(tmp_pdb_path)
                
    except Exception as e:
        logger.error(f"Error converting PDB to SDF: {e}")
        return None


def pdb_file_to_sdf_file(pdb_path: str, sdf_path: str, 
                         remove_waters: bool = True, remove_ions: bool = True) -> bool:
    """
    Convert a PDB file to SDF file.
    
    Args:
        pdb_path: Path to input PDB file
        sdf_path: Path to output SDF file
        remove_waters: Remove water molecules
        remove_ions: Remove common ions
        
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()
        
        sdf_content = pdb_to_sdf(pdb_content, remove_waters, remove_ions)
        
        if sdf_content:
            with open(sdf_path, 'w') as f:
                f.write(sdf_content)
            return True
        return False
        
    except Exception as e:
        logger.error(f"Error converting PDB file to SDF: {e}")
        return False


def ligand_atomgroup_to_sdf(ligand_atomgroup, res_name: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Convert a ProDy AtomGroup ligand to SDF format with proper bond orders.
    
    Args:
        ligand_atomgroup: ProDy AtomGroup containing the ligand
        res_name: Residue name of the ligand
        
    Returns:
        Tuple of (SDF content, SMILES string) or (None, None) if conversion failed
    """
    try:
        # Process ligand with proper bond order assignment
        mols = process_ligand_stereochemistry(ligand_atomgroup, res_name)
        
        if not mols:
            logger.error(f"Failed to process ligand {res_name}")
            return None, None
        
        # Use the first conformer
        mol = mols[0]
        
        # Generate SDF content
        sdf_content = Chem.MolToMolBlock(mol)
        
        # Generate SMILES
        smiles = Chem.MolToSmiles(mol)
        
        return sdf_content, smiles
        
    except Exception as e:
        logger.error(f"Error converting ligand AtomGroup to SDF: {e}")
        return None, None


def extract_smiles_from_pdb(pdb_content: str, res_name: Optional[str] = None,
                           remove_waters: bool = True, remove_ions: bool = True) -> Optional[str]:
    """
    Extract SMILES string from PDB content.
    
    Args:
        pdb_content: PDB file content as string
        remove_waters: Remove water molecules
        remove_ions: Remove common ions
        
    Returns:
        SMILES string, or None if extraction failed
    """
    try:
        # Write PDB content to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
            tmp_pdb.write(pdb_content)
            tmp_pdb_path = tmp_pdb.name
        
        try:
            # Read PDB file with RDKit
            mol = Chem.MolFromPDBFile(tmp_pdb_path, removeHs=False, sanitize=True)
            
            if mol is None:
                # Try without sanitization
                mol = Chem.MolFromPDBFile(tmp_pdb_path, removeHs=False, sanitize=False)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                    except:
                        pass
            
            if mol is None:
                return None
            
            # Remove waters and ions if requested
            if remove_waters or remove_ions:
                editable_mol = Chem.EditableMol(mol)
                atoms_to_remove = []
                
                for atom in mol.GetAtoms():
                    residue_info = atom.GetPDBResidueInfo()
                    if residue_info:
                        res_name = residue_info.GetResidueName().strip()
                        if remove_waters and res_name in ['HOH', 'WAT', 'H2O']:
                            atoms_to_remove.append(atom.GetIdx())
                        elif remove_ions and res_name in ['NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE']:
                            atoms_to_remove.append(atom.GetIdx())
                
                # Remove atoms in reverse order
                for idx in sorted(atoms_to_remove, reverse=True):
                    editable_mol.RemoveAtom(idx)
                
                mol = editable_mol.GetMol()
            
            # Generate SMILES
            smiles = Chem.MolToSmiles(mol)
            return smiles
            
        finally:
            # Clean up temp file
            if os.path.exists(tmp_pdb_path):
                os.unlink(tmp_pdb_path)
                
    except Exception as e:
        logger.error(f"Error extracting SMILES from PDB: {e}")
        return None


def extract_smiles_from_sdf(sdf_content: str) -> Optional[str]:
    """
    Extract SMILES string from SDF content.
    
    Args:
        sdf_content: SDF file content as string
        
    Returns:
        SMILES string, or None if extraction failed
    """
    try:
        mol = Chem.MolFromMolBlock(sdf_content)
        if mol:
            return Chem.MolToSmiles(mol)
        return None
    except Exception as e:
        logger.error(f"Error extracting SMILES from SDF: {e}")
        return None





def calculate_ligand_properties(smiles: str) -> dict:
    """
    Calculate molecular properties for a ligand from its SMILES.
    
    Args:
        smiles: SMILES string
        
    Returns:
        Dictionary of molecular properties
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        properties = {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "tpsa": Descriptors.TPSA(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": Descriptors.RingCount(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
            "num_heteroatoms": Descriptors.NumHeteroatoms(mol),
            "formal_charge": Chem.rdmolops.GetFormalCharge(mol)
        }
        
        return properties
        
    except Exception as e:
        logger.error(f"Error calculating ligand properties: {e}")
        return {}


def split_multi_molecule_sdf(sdf_content: str) -> list:
    """
    Split a multi-molecule SDF into individual molecules.
    
    Args:
        sdf_content: SDF file content potentially containing multiple molecules
        
    Returns:
        List of individual SDF contents
    """
    molecules = []
    try:
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_content)
        
        for mol in suppl:
            if mol:
                molecules.append(Chem.MolToMolBlock(mol))
    except Exception as e:
        logger.error(f"Error splitting multi-molecule SDF: {e}")
    
    return molecules


def AssignBondOrdersFromTemplate(refmol: Chem.Mol, mol: Chem.Mol) -> list[Chem.Mol]:
    """
    Assign bond orders from a template to an RDKit molecule, handling all possible matches
    :param refmol: template molecule
    :param mol: RDKit molecule
    :return: List of RDKit molecules with bond orders assigned for each possible match
    """
    
    refmol2 = Chem.Mol(refmol)
    mol2 = Chem.Mol(mol)
    result_molecules = []
    
    # do the molecules match already?
    matching = mol2.GetSubstructMatch(refmol2)
    if not matching:  # no, they don't match
        # check if bonds of mol are SINGLE
        for b in mol2.GetBonds():
            if b.GetBondType() != Chem.BondType.SINGLE:
                b.SetBondType(Chem.BondType.SINGLE)
                b.SetIsAromatic(False)
        # set the bonds of mol to SINGLE
        for b in refmol2.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        # set atom charges to zero;
        for a in refmol2.GetAtoms():
            a.SetFormalCharge(0)
        for a in mol2.GetAtoms():
            a.SetFormalCharge(0)

        matchings = mol2.GetSubstructMatches(refmol2, uniquify=False)
        # do the molecules match now?
        if matchings:
            for matching in matchings:
                # Create a new copy of the molecule for each match
                current_mol = Chem.Mol(mol2)
                # apply matching: set bond properties
                for b in refmol.GetBonds():
                    atom1 = matching[b.GetBeginAtomIdx()]
                    atom2 = matching[b.GetEndAtomIdx()]
                    b2 = current_mol.GetBondBetweenAtoms(atom1, atom2)
                    b2.SetBondType(b.GetBondType())
                    b2.SetIsAromatic(b.GetIsAromatic())
                # apply matching: set atom properties
                for a in refmol.GetAtoms():
                    a2 = current_mol.GetAtomWithIdx(matching[a.GetIdx()])
                    a2.SetHybridization(a.GetHybridization())
                    a2.SetIsAromatic(a.GetIsAromatic())
                    a2.SetNumExplicitHs(a.GetNumExplicitHs())
                    a2.SetFormalCharge(a.GetFormalCharge())
                Chem.SanitizeMol(current_mol)
                if hasattr(current_mol, '__sssAtoms'):
                    current_mol.__sssAtoms = None  # we don't want all bonds highlighted
                result_molecules.append(current_mol)
        else:
            raise ValueError("No matching found")
    else:
        # If there's an exact match, just return the original molecule
        result_molecules.append(mol2)
        
    return result_molecules


def process_ligand_stereochemistry(ligand, res_name):
    """
    Add bond orders to a pdb ligand (ProDy AtomGroup)
    1. Select the ligand component with name "res_name"
    2. Get the corresponding SMILES from pypdb
    3. Create a template molecule from the SMILES in step 2
    4. Write the PDB file to a stream
    5. Read the stream into an RDKit molecule
    6. Assign the bond orders from the template from step 3
    :param ligand: ligand as generated by prody
    :param res_name: residue name of ligand to extract
    :return: list of molecules with bond orders assigned
    """
    output = StringIO()
    sub_mol = ligand.select(f"resname {res_name}")
    chem_desc = pypdb.describe_chemical(f"{res_name}")
    sub_smiles = chem_desc["rcsb_chem_comp_descriptor"]["smiles"]
    template = AllChem.MolFromSmiles(sub_smiles)
    writePDBStream(output, sub_mol)
    pdb_string = output.getvalue()
    rd_mol = AllChem.MolFromPDBBlock(pdb_string)
    new_mols = AssignBondOrdersFromTemplate(template, rd_mol)
    for idx, mol in enumerate(new_mols):
        # Add hydrogens and generate coordinates only for the new H atoms.
        # Setting addCoords=True preserves all existing (heavy‑atom) coordinates.
        mol_with_h = Chem.AddHs(mol, addCoords=True)
        new_mols[idx] = mol_with_h
        
    return new_mols


def write_protein_pdbqt(protein_path: str, outfile: str, allow_bad_res: bool = True) -> str:
    """
    Convert protein PDB to PDBQT format for docking using Meeko.
    
    Args:
        protein_path: Path to protein PDB file
        outfile: Path for output PDBQT file
        allow_bad_res: Whether to allow non-standard residues
    
    Returns:
        Path to the created PDBQT file
    """
    import subprocess
    
    cmd = ["mk_prepare_receptor.py", "-i", protein_path, "-o", outfile]
    if allow_bad_res:
        cmd.append("--allow_bad_res")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Failed to prepare protein PDBQT: {result.stderr}")
    
    if not os.path.exists(outfile):
        raise FileNotFoundError(f"PDBQT file was not created: {outfile}")
    
    logger.info(f"Created protein PDBQT: {outfile}")
    return outfile


def write_ligand_pdbqt(ligand_path: str, outfile: str, input_format: str = "sdf") -> str:
    """
    Convert ligand to PDBQT format for docking using Meeko.
    
    Args:
        ligand_path: Path to ligand file (SDF, MOL2, or PDB)
        outfile: Path for output PDBQT file
        input_format: Input file format (sdf, mol2, pdb)
    
    Returns:
        Path to the created PDBQT file
    """
    import subprocess
    
    cmd = ["mk_prepare_ligand.py", "-i", ligand_path, "-o", outfile]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Failed to prepare ligand PDBQT: {result.stderr}")
    
    if not os.path.exists(outfile):
        raise FileNotFoundError(f"PDBQT file was not created: {outfile}")
    
    logger.info(f"Created ligand PDBQT: {outfile}")
    return outfile


def add_missing_protein_atoms(pdb_path: str, output_path: Optional[str] = None, 
                             add_hydrogens: bool = True, ph: float = 7.0) -> str:
    """
    Use PDBFixer to add missing atoms and replace non-standard residues.
    
    Args:
        pdb_path: Path to input PDB file
        output_path: Path for output PDB file (if None, overwrites input)
        add_hydrogens: Whether to add missing hydrogens
        ph: pH for hydrogen addition (default 7.0)
    
    Returns:
        Path to the fixed PDB file
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    
    fixer = PDBFixer(pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    
    if add_hydrogens:
        # add hydrogens to the protein at specified pH
        fixer.addMissingHydrogens(ph)
    
    # Determine output path
    if output_path is None:
        output_path = pdb_path
    
    with open(output_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    
    return output_path

def vina_score(protein, ligand, center: list[float], box_size: list[float]) -> float:
    """
    Score a ligand using Vina
    :param protein: path to protein pdbqt file
    :param ligand: path to ligand pdbqt file
    :param center: center coordinates for docking box
    :param box_size: size of docking box
    :return: score
    """
    import os
    from vina import Vina

    # Verify files exist
    if not os.path.exists(protein):
        raise FileNotFoundError(f"Protein file not found: {protein}")
    if not os.path.exists(ligand):
        raise FileNotFoundError(f"Ligand file not found: {ligand}")

    v = Vina(sf_name='vina')
    v.set_receptor(protein)
    v.set_ligand_from_file(ligand)
    
    v.compute_vina_maps(center=center, box_size=box_size)

    # Score the current pose
    energy = v.score()
    return energy

