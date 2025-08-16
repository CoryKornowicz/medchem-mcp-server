import base64
import io
from PIL import Image
import numpy as np
import prody
from scipy.spatial import cKDTree
from prody import Bond

# DEFINE TOOLS
def pil_image_to_base64(img, format="PNG"):
    """
    Converts a PIL Image object to a base64 string.

    Args:
        img: PIL Image object.
        format: Image format for saving to buffer (e.g., "PNG", "JPEG").

    Returns:
        A base64 encoded string representation of the image.
    """
    buffered = io.BytesIO()
    img.save(buffered, format=format)
    img_byte = buffered.getvalue()
    return base64.b64encode(img_byte).decode()

def img_byte_to_base64(img_byte):
    """
    Converts a byte array to a base64 string.

    Args:
        img_byte: Byte array of the image.

    Returns:
        A base64 encoded string representation of the image.
    """
    return base64.b64encode(img_byte).decode()


def fastInferBonds_kdtree_vectorized(structure: prody.atomic.atomgroup.AtomGroup, max_bond=1.6, min_bond=0, set_bonds=True) -> list:
    """
    Even faster vectorized implementation using KDTree and numpy operations.
    """
    
    coords = structure.getCoords()
    n_atoms = structure._n_atoms
    indices = np.arange(structure._n_atoms)
    
    if n_atoms < 2:
        return []
    
    # Build KDTree
    kdtree = cKDTree(coords)
    
    # Find all pairs within max_bond distance
    pairs = kdtree.query_pairs(max_bond, output_type='ndarray')
    
    if len(pairs) == 0:
        return []
    
    # Calculate distances for all pairs at once
    coord_pairs = coords[pairs]
    distances = np.linalg.norm(coord_pairs[:, 0] - coord_pairs[:, 1], axis=1)
    
    # Filter by min_bond distance
    valid_pairs = pairs[distances > min_bond]
    
    # Convert to bond format using original indices
    bonds = [[indices[pair[0]], indices[pair[1]]] for pair in valid_pairs]
    
    if set_bonds:
        structure.setBonds(bonds)

    # return bonds
    acsi = structure._acsi
    return np.array([Bond(structure, bond, acsi) for bond in bonds])

# Determining functions 

def determine_non_standard_amino_acids(structure):
    """
    Identify non-standard amino acids in a structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The protein structure to analyze
        
    Returns:
    --------
    set: Names of non-standard amino acids found in the structure
    """
    non_standard_aa = set()
    
    for res in structure.iterResidues():
        if _looks_like_non_standard_amino_acid(res):
            non_standard_aa.add(res.getResname())

    if len(non_standard_aa) == 0:
        return None
    # extract from structure non-standard amino acids
    non_standard_aa_atoms = structure.select(f'resname {" ".join(non_standard_aa)}')
    return non_standard_aa_atoms

def _looks_like_non_standard_amino_acid(res) -> bool:
    """Return True if the residue is HETERO *and* has a peptide backbone."""
    # Any atom flagged HETERO → residue is hetero
    hetero = any(atom.getFlag("hetero") for atom in res)          # Atom objects do have .isHet()
    # Backbone check: N, CA, C all present
    backbone = all(a in res.getNames() for a in ("N", "CA", "C"))
    return hetero and backbone

def determine_ligands_in_system(structure):
    """
    Comprehensive function to identify ligands in a protein system.
    Returns indices of atoms that are ligands (small molecules, disconnected amino acids, etc.)
    
    This function properly detects:
    1. Small molecules and cofactors (non-standard residues without protein backbone)
    2. Isolated amino acids (like disconnected TYR, etc.)
    3. Terminal residues with OXT that are disconnected
    4. Small disconnected nucleic acids
    5. Non-standard residues that are truly disconnected from the main protein structure
    
    This function SHOULD NOT detect: 
    1. Waters (HOH, WAT, etc.)
    2. Ions (like Na+, Cl-, etc.)
    
    The function is conservative about non-standard amino acids and uses bond connectivity
    analysis to determine if residues are truly disconnected from the main protein structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The protein structure to analyze
        
    Returns:
    --------
    list: Indices of atoms that are ligands
    """
    
    # Standard amino acid residue names
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    # Standard nucleic acid residue names
    standard_na = {'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC', 'DU'}
    
    # Common solvent molecules to exclude
    solvents = {'HOH', 'WAT', 'TIP', 'SOL', 'H2O'}
    
    # Standard ions to exclude
    standard_ions = {'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'CU', 'MN'}
    
    ligand_indices = []
    
    # First, build a chain connectivity map
    chain_connectivity = {}
    for res in structure.iterResidues():
        chain = res.getChid()
        resnum = res.getResnum()
        resname = res.getResname()
        
        if chain not in chain_connectivity:
            chain_connectivity[chain] = {}
        chain_connectivity[chain][resnum] = resname
    
    for res in structure.iterResidues():
        resname = res.getResname()
        chain = res.getChid()
        resnum = res.getResnum()
        
        # Skip solvent molecules
        if resname in solvents or resname in standard_ions:
            continue
            
        # Check if it's a hetero residue (non-standard)
        is_hetero = any(atom.getFlag("hetero") for atom in res)
        
        # Case 1: Non-standard residues - be more conservative
        if resname not in standard_aa and resname not in standard_na:
            # Check if it looks like a modified amino acid with proper backbone
            has_aa_backbone = all(name in res.getNames() for name in ['N', 'CA', 'C'])
            
            # If it has AA backbone, check if it's truly disconnected
            if has_aa_backbone:
                # Only consider it a ligand if it's both small AND disconnected
                # This prevents modified amino acids in the main chain from being filtered
                if len(res) < 25 and _is_disconnected_from_main_chain(res, structure):
                    print(f"Identified ligand-like residue: {resname} at {chain}:{resnum}")
                    print(f"Disconnected from main chain, {_is_disconnected_from_main_chain(res, structure)}.")
                    ligand_indices.extend(res.getIndices())
                # # Also check if it's isolated (no sequential neighbors)
                # elif _is_residue_isolated(res, chain_connectivity):
                #     ligand_indices.extend(res.getIndices())
            else:
                # No standard backbone - likely a true ligand (cofactor, small molecule, etc.)
                # But still check if it's not just a large polymer connected to the protein
                # TODO: Test if this allows ACE, NME, NH2, and other connected fragments | len(res) < 50 or
                if _is_disconnected_from_main_chain(res, structure):
                    ligand_indices.extend(res.getIndices())
        
        # Case 2: Standard amino acids that might be ligands (isolated/disconnected)
        elif resname in standard_aa:
            # Check for isolation (key improvement!)
            is_isolated = _is_residue_isolated(res, chain_connectivity)
            
            # Check if it has terminal markers
            has_oxt = 'OXT' in res.getNames()
            
            # If it's hetero-flagged, isolated, or has terminal markers, it's likely a ligand
            if is_hetero or (is_isolated and has_oxt) or (is_isolated and len(res) < 20):
                ligand_indices.extend(res.getIndices())
        
        # Case 3: Standard nucleic acids that are disconnected
        elif resname in standard_na:
            if is_hetero or _is_disconnected_from_main_chain(res, structure):
                ligand_indices.extend(res.getIndices())
    
    if len(ligand_indices) == 0:
        return None
    
    return ligand_indices

def _is_residue_isolated(residue, chain_connectivity):
    """
    Check if a residue is isolated (no sequential neighbors) in its chain.
    
    Parameters:
    -----------
    residue : prody.Residue
        The residue to check
    chain_connectivity : dict
        Dictionary mapping chain -> {resnum: resname}
        
    Returns:
    --------
    bool: True if residue has no sequential neighbors
    """
    chain = residue.getChid()
    resnum = residue.getResnum()
    
    if chain not in chain_connectivity:
        return True
    
    chain_residues = chain_connectivity[chain]
    
    # Check if previous (resnum-1) and next (resnum+1) residues exist
    has_prev = (resnum - 1) in chain_residues
    has_next = (resnum + 1) in chain_residues
    
    # If both neighbors are missing, it's isolated
    return not has_prev and not has_next

def _is_disconnected_from_main_chain(residue, structure, max_bond_distance=1.8):
    """
    Check if a residue is disconnected from the main protein chain using bond connectivity.
    
    Parameters:
    -----------
    residue : prody.Residue
        The residue to check
    structure : prody.AtomGroup  
        The full structure
    max_bond_distance : float
        Maximum distance for considering atoms bonded (default 1.8 Å)
        
    Returns:
    --------
    bool: True if residue appears disconnected from main chain
    """
    try:
        # Get all atoms in the structure
        all_atoms = structure.select('all')
        if all_atoms is None:
            return True

        # Build bonds for the entire structure if not already present
        bonds = structure.getBonds()
        if bonds is None:
            try:
                # bonds = structure.inferBonds()
                bonds = fastInferBonds_kdtree_vectorized(structure, set_bonds=True)
            except:
                print("Failed to build bond list for structure. Falling back to distance-based method.")
                # If bond building fails, fall back to distance-based method
                return _distance_based_connectivity_check(residue, structure)

        # Get atoms from this residue
        res_atoms = residue.select('all')
        if res_atoms is None:
            return True

        # Get all other residues (excluding current residue)
        other_residues = []
        for other_res in structure.iterResidues():
            if (other_res.getResnum() != residue.getResnum()):
                other_residues.append(other_res)

        if not other_residues:
            return True

        # Check if any atom in the residue is bonded to atoms in other residues
        res_indices = set(res_atoms.getIndices())
        for bond in bonds:
            atom1_idx, atom2_idx = bond
            
            # Check if one atom is in our residue and the other is not
            atom1_in_res = atom1_idx.getIndex() in res_indices
            atom2_in_res = atom2_idx.getIndex() in res_indices

            if atom1_in_res and not atom2_in_res:
                # Found a bond from residue to outside - it's connected
                return False
            elif atom2_in_res and not atom1_in_res:
                # Found a bond from residue to outside - it's connected
                return False

        # No bonds found to other residues - it's disconnected
        return True
       
    except Exception as e:
        print(f"Error during bond connectivity check: {e}")
        # If bond analysis fails, fall back to distance-based method
        return _distance_based_connectivity_check(residue, structure)

def _distance_based_connectivity_check(residue, structure, distance_cutoff=2.0):
    """
    Fallback distance-based connectivity check.
    """
    try:
        from prody import calcDistance
        
        # Get protein atoms (standard amino acids) excluding this residue
        protein_selection = f'protein and not (chain {residue.getChid()} and resnum {residue.getResnum()})'
        protein_atoms = structure.select(protein_selection)
        
        if protein_atoms is None:
            return True
            
        # Get atoms from this residue
        res_atoms = residue.select('all')
        
        # Sample a few atoms from the residue for efficiency
        sample_atoms = res_atoms[:min(5, len(res_atoms))]
        
        for atom in sample_atoms:
            # Find closest protein atom
            distances = calcDistance(atom, protein_atoms)
            min_distance = min(distances) if len(distances) > 0 else float('inf')
            
            # If any atom is close to protein, consider it connected
            if min_distance < distance_cutoff:
                return False
                
        return True
        
    except:
        # If there's any error, assume it's disconnected
        return True

def get_ions(structure):
    """
    Extract ions from a structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The structure to analyze
        
    Returns:
    --------
    prody.AtomGroup: Ions in the structure
    """
    
    # Standard ions to exclude
    standard_ions = {'NA', 'CL', 'MG', 'CA', 'K', 'ZN', 'FE', 'CU', 'MN'}

    ion_selection = structure.select(f'resname {" ".join(standard_ions)}')

    if ion_selection is None:
        return None

    return structure[ion_selection.getIndices()].toAtomGroup()

def get_nucleic_acids(structure):
    """Extract nucleic acids from a structure.
    Parameters:
    
    -----------
    structure : prody.AtomGroup
        The structure to analyze
    Returns:
    --------
    prody.AtomGroup: Nucleic acids in the structure
    """
    # Select nucleic acids (DNA/RNA) based on standard residue names
    nucleic_selection = structure.select('nucleic')

    if nucleic_selection is None:
        return None

    return structure[nucleic_selection.getIndices()].toAtomGroup()

def get_ligands(structure):
    """
    Extract only the ligands from a structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The structure to analyze
        
    Returns:
    --------
    prody.AtomGroup: Only the ligand atoms
    """
    ligand_indices = determine_ligands_in_system(structure)

    if ligand_indices is None:
        return None

    return structure[ligand_indices].toAtomGroup()

def get_non_standard_amino_acids(structure):
    """
    Get non-standard amino acids from a structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The structure to analyze
        
    Returns:
    --------
    prody.AtomGroup: Non-standard amino acids found in the structure
    """    
    nn_aa_indices = determine_non_standard_amino_acids(structure)

    if nn_aa_indices is None:
        return None

    return structure[nn_aa_indices.getIndices()].toAtomGroup()

def get_waters(structure):
    """
    Extract water molecules from a structure.
    
    Parameters:
    -----------
    structure : prody.AtomGroup
        The structure to analyze
        
    Returns:
    --------
    prody.AtomGroup: Water molecules in the structure
    """
    
    water_sel = structure.select('water')
    
    if water_sel is None:
        return None

    return water_sel.toAtomGroup()

def get_protein(structure):
    """
    Extract protein molecules from a structure.
    This includes retaining non-standard amino acids that are part of the protein structure.

    Parameters:
    -----------
    structure : prody.AtomGroup
        The structure to analyze

    Returns:
    --------
    prody.AtomGroup: Protein molecules in the structure
    """
    # Select all protein atoms, including non-standard amino acids
    protein_selection = 'protein'
    
    nn_res = get_non_standard_amino_acids(structure)
    if nn_res is not None:
        protein_selection = 'protein or resname ' + ' or resname '.join(f'"{res}"' for res in nn_res)
    
    ligands = get_ligands(structure)
    if ligands is not None:
        protein_selection = f"({protein_selection}) and not resname {" ".join(str(res) for res in set(ligands.getResnames()))}"

    protein_group = structure.select(protein_selection)
    return protein_group.toAtomGroup()

#  EXAMPLES --------- 
# structure = parsePDB(test_file)
# if structure is None:
#     print(f"Failed to parse {test_file}")
#     continue
# print(f"Processing {test_file}...")
# # get proteins, ligands, waters, ions, nucleic acids, non-standard amino acids
# prot = get_protein(structure)
# ligands = get_ligands(structure)
# waters = get_waters(structure)
# ions = get_ions(structure)
# nucleic_acids = get_nucleic_acids(structure)
# non_standard_aa = get_non_standard_amino_acids(structure)

# print(f"Protein: {len(prot)} atoms")
# if ligands is not None:
#     print(f"Ligands: {len(ligands)} atoms")
#     print(set(ligands.getResnames()))
# if waters is not None:
#     print(f"Waters: {len(waters)} atoms")
#     print(set(waters.getResnames()))
# if ions is not None:
#     print(f"Ions: {len(ions)} atoms")
#     print(set(ions.getResnames()))
# if nucleic_acids is not None:
#     print(f"Nucleic Acids: {len(nucleic_acids)} atoms")
#     print(set(nucleic_acids.getResnames()))
# if non_standard_aa is not None:
#     print(f"Non-standard Amino Acids: {len(non_standard_aa)} atoms")
#     print(set(non_standard_aa.getResnames()))
#     # check if the non-standard amino acid is disconnected from the main chain
#     for res in non_standard_aa.iterResidues():
#         disconnected = _is_disconnected_from_main_chain(res, structure)
#         if disconnected:
#             print(f"Non-standard amino acid {res.getResname()} at {res.getChid()}:{res.getResnum()} is disconnected from the main chain.")
#         else:
#             print(f"Non-standard amino acid {res.getResname()} at {res.getChid()}:{res.getResnum()} is connected to the main chain.")

# OUTPUT --------- 
# Processing data/1hsg.pdb...
# Protein: 1514 atoms
# Ligands: 45 atoms
# {'MK1'}
# Waters: 127 atoms
# {'HOH'}
# Processing data/1ubq.pdb...
# Protein: 602 atoms
# Waters: 58 atoms
# {'HOH'}
# Processing data/4oud.pdb...
# Protein: 3100 atoms
# Ligands: 13 atoms
# {'TYR'}
# Waters: 37 atoms
# {'HOH'}
# Non-standard Amino Acids: 17 atoms
# {'BIF'}
# Non-standard amino acid BIF at A:303 is connected to the main chain.
# Processing data/5kgn.pdb...
# Protein: 8147 atoms
# Ligands: 12 atoms
# {'GOL'}
# Waters: 735 atoms
# {'HOH'}
# Ions: 8 atoms
# {'CL', 'ZN', 'MN', 'MG'}
# Non-standard Amino Acids: 24 atoms
# {'DTY'}
# Non-standard amino acid DTY at C:1 is connected to the main chain.
# Non-standard amino acid DTY at D:1 is connected to the main chain.
# Processing data/alanine-dipeptide.pdb...
# Protein: 10 atoms

