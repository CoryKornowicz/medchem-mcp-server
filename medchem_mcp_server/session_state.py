"""
Session state management for the medchem MCP server.

This module provides functionality to store and manage ephemeral data
per session, such as molecule collections, protein structures, and files
that don't need to be returned to the chat context.
"""

import os
import uuid
import tempfile
import shutil
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field
from datetime import datetime
import threading
from pathlib import Path

# Lazy import for protein analysis utilities to avoid dependency issues
UTILITIES_AVAILABLE = False
_utilities_imported = False

def _import_utilities():
    """Lazy import utilities when needed"""
    global UTILITIES_AVAILABLE, _utilities_imported
    global get_protein, get_ligands, get_waters, get_ions, get_nucleic_acids, get_non_standard_amino_acids
    
    if _utilities_imported:
        return UTILITIES_AVAILABLE
    
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
        _utilities_imported = True
    except ImportError:
        UTILITIES_AVAILABLE = False
        _utilities_imported = True
    
    return UTILITIES_AVAILABLE

@dataclass
class MoleculeEntry:
    """Represents a stored molecule with metadata"""
    smiles: str
    name: Optional[str] = None
    properties: Dict[str, Any] = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)
    tags: List[str] = field(default_factory=list)

@dataclass
class ProteinEntry:
    """Represents a stored protein structure with metadata"""
    pdb_id: str
    name: Optional[str] = None
    pdb_content: Optional[str] = None  # Raw PDB file content
    file_path: Optional[str] = None  # Path to saved PDB file
    chain_count: Optional[int] = None
    residue_count: Optional[int] = None
    properties: Dict[str, Any] = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)
    tags: List[str] = field(default_factory=list)

@dataclass
class FileEntry:
    """Represents a stored file with metadata"""
    file_id: str
    original_name: str
    file_path: str  # Absolute path to file in tmp directory
    file_type: str  # e.g., 'pdb', 'sdf', 'mol2', 'csv', etc.
    size_bytes: int
    properties: Dict[str, Any] = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)
    tags: List[str] = field(default_factory=list)

@dataclass 
class SessionData:
    """Container for session-specific data"""
    session_id: str
    molecules: Dict[str, MoleculeEntry] = field(default_factory=dict)
    proteins: Dict[str, ProteinEntry] = field(default_factory=dict)
    files: Dict[str, FileEntry] = field(default_factory=dict)
    collections: Dict[str, List[str]] = field(default_factory=dict)  # collection_name -> entity_ids
    tmp_dir: Optional[str] = None  # Session-specific temp directory
    metadata: Dict[str, Any] = field(default_factory=dict)
    created_at: datetime = field(default_factory=datetime.now)
    last_accessed: datetime = field(default_factory=datetime.now)

class SessionStateManager:
    """Manages ephemeral state per session"""
    
    def __init__(self):
        self._sessions: Dict[str, SessionData] = {}
        self._lock = threading.Lock()
        # Create base temp directory for all sessions
        self._base_tmp_dir = tempfile.mkdtemp(prefix="medchem_mcp_")
    
    def __del__(self):
        """Cleanup all temp directories on shutdown"""
        try:
            if hasattr(self, '_base_tmp_dir') and os.path.exists(self._base_tmp_dir):
                shutil.rmtree(self._base_tmp_dir)
        except Exception:
            pass  # Best effort cleanup
    
    def get_or_create_session(self, session_id: Optional[str] = None) -> str:
        """Get existing session or create new one"""
        if session_id is None:
            session_id = str(uuid.uuid4())
        
        with self._lock:
            if session_id not in self._sessions:
                # Create session-specific temp directory
                session_tmp_dir = os.path.join(self._base_tmp_dir, session_id)
                os.makedirs(session_tmp_dir, exist_ok=True)
                
                self._sessions[session_id] = SessionData(
                    session_id=session_id,
                    tmp_dir=session_tmp_dir
                )
            else:
                self._sessions[session_id].last_accessed = datetime.now()
        
        return session_id
    
    def get_session(self, session_id: str) -> Optional[SessionData]:
        """Get session data"""
        with self._lock:
            session = self._sessions.get(session_id)
            if session:
                session.last_accessed = datetime.now()
            return session
    
    # Molecule management methods
    def add_molecule(self, session_id: str, smiles: str, name: Optional[str] = None, 
                    properties: Optional[Dict[str, Any]] = None, tags: Optional[List[str]] = None) -> str:
        """Add molecule to session and return molecule ID"""
        molecule_id = str(uuid.uuid4())
        molecule = MoleculeEntry(
            smiles=smiles,
            name=name,
            properties=properties or {},
            tags=tags or []
        )
        
        session = self.get_session(session_id)
        if session:
            session.molecules[molecule_id] = molecule
        
        return molecule_id
    
    def get_molecule(self, session_id: str, molecule_id: str) -> Optional[MoleculeEntry]:
        """Get molecule from session"""
        session = self.get_session(session_id)
        return session.molecules.get(molecule_id) if session else None
    
    def list_molecules(self, session_id: str, tags: Optional[List[str]] = None) -> Dict[str, MoleculeEntry]:
        """List molecules in session, optionally filtered by tags"""
        session = self.get_session(session_id)
        if not session:
            return {}
        
        if tags is None:
            return session.molecules.copy()
        
        filtered = {}
        for mol_id, molecule in session.molecules.items():
            if any(tag in molecule.tags for tag in tags):
                filtered[mol_id] = molecule
        
        return filtered
    
    # Protein management methods
    def add_protein(self, session_id: str, pdb_id: str, name: Optional[str] = None,
                   pdb_content: Optional[str] = None, file_path: Optional[str] = None,
                   properties: Optional[Dict[str, Any]] = None, tags: Optional[List[str]] = None) -> str:
        """Add protein to session and return protein ID"""
        protein_id = str(uuid.uuid4())
        
        # Parse basic protein info if content provided using ProDy and utilities
        chain_count = None
        residue_count = None
        atom_count = None
        chain_info = {}
        component_info = {}
        
        if pdb_content or file_path:
            try:
                import prody
                import tempfile
                import os
                
                # Determine path to PDB file
                if file_path and os.path.exists(file_path):
                    parse_path = file_path
                    need_cleanup = False
                elif pdb_content:
                    # Write content to temp file for ProDy parsing
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
                        tmp.write(pdb_content)
                        parse_path = tmp.name
                    need_cleanup = True
                else:
                    parse_path = None
                    need_cleanup = False
                
                if parse_path:
                    try:
                        # Parse with ProDy
                        structure = prody.parsePDB(parse_path)
                        
                        if structure:
                            # Use utility functions if available for accurate segmentation
                            if _import_utilities():
                                # Get session temp directory for saving extracted components
                                session = self.get_session(session_id)
                                if session and session.tmp_dir:
                                    # Get protein component specifically
                                    protein_atoms = get_protein(structure)
                                    if protein_atoms:
                                        # Get chain information from protein component only
                                        chains = set(protein_atoms.getChids())
                                        chain_count = len(chains)
                                        
                                        # Count residues per chain in protein component
                                        for chain_id in chains:
                                            chain_sel = protein_atoms.select(f'chain {chain_id}').toAtomGroup()
                                            if chain_sel:
                                                chain_residues = set()
                                                for res in chain_sel.iterResidues():
                                                    chain_residues.add((res.getChid(), res.getResnum()))
                                                chain_info[chain_id] = len(chain_residues)
                                        
                                        # Total residue count in protein
                                        all_residues = set()
                                        for res in protein_atoms.iterResidues():
                                            all_residues.add((res.getChid(), res.getResnum()))
                                        residue_count = len(all_residues)
                                        
                                        component_info['protein_atoms'] = len(protein_atoms)
                                        
                                        # Save extracted protein to separate file
                                        protein_only_path = os.path.join(session.tmp_dir, f"{pdb_id}_protein.pdb")
                                        prody.writePDB(protein_only_path, protein_atoms)
                                        component_info['protein_file'] = protein_only_path
                                    
                                    # Get ligands and save separately
                                    ligands = get_ligands(structure)
                                    if ligands:
                                        component_info['ligand_atoms'] = len(ligands)
                                        component_info['ligand_names'] = list(set(ligands.getResnames()))
                                        
                                        # Save each ligand type to separate file
                                        ligand_files = {}
                                        for ligand_name in component_info['ligand_names']:
                                            ligand_sel = ligands.select(f'resname {ligand_name}').toAtomGroup()
                                            if ligand_sel:
                                                ligand_path = os.path.join(session.tmp_dir, f"{pdb_id}_{ligand_name}.pdb")
                                                prody.writePDB(ligand_path, ligand_sel)
                                                ligand_files[ligand_name] = ligand_path
                                        
                                        # Also save all ligands together
                                        all_ligands_path = os.path.join(session.tmp_dir, f"{pdb_id}_ligands.pdb")
                                        prody.writePDB(all_ligands_path, ligands)
                                        component_info['ligand_files'] = ligand_files
                                        component_info['all_ligands_file'] = all_ligands_path
                                    
                                    # Get waters
                                    waters = get_waters(structure)
                                    if waters:
                                        water_residues = set()
                                        for res in waters.iterResidues():
                                            water_residues.add((res.getChid(), res.getResnum()))
                                        component_info['water_molecules'] = len(water_residues)
                                        
                                        # Optionally save waters to file
                                        waters_path = os.path.join(session.tmp_dir, f"{pdb_id}_waters.pdb")
                                        prody.writePDB(waters_path, waters)
                                        component_info['waters_file'] = waters_path
                                    
                                    # Get ions
                                    ions = get_ions(structure)
                                    if ions:
                                        component_info['ion_atoms'] = len(ions)
                                        component_info['ion_types'] = list(set(ions.getResnames()))
                                        
                                        # Save ions to file
                                        ions_path = os.path.join(session.tmp_dir, f"{pdb_id}_ions.pdb")
                                        prody.writePDB(ions_path, ions)
                                        component_info['ions_file'] = ions_path
                                    
                                    # Get non-standard amino acids
                                    non_standard_aa = get_non_standard_amino_acids(structure)
                                    if non_standard_aa:
                                        component_info['non_standard_aa_atoms'] = len(non_standard_aa)
                                        component_info['non_standard_aa_types'] = list(set(non_standard_aa.getResnames()))
                                
                            else:
                                # Fallback to basic ProDy analysis
                                chains = set(structure.getChids())
                                chain_count = len(chains)
                                
                                for chain_id in chains:
                                    chain_sel = structure.select(f'chain {chain_id}').toAtomGroup()
                                    if chain_sel:
                                        chain_residues = set()
                                        for res in chain_sel.iterResidues():
                                            chain_residues.add((res.getChid(), res.getResnum()))
                                        chain_info[chain_id] = len(chain_residues)
                                
                                all_residues = set()
                                for res in structure.iterResidues():
                                    all_residues.add((res.getChid(), res.getResnum()))
                                residue_count = len(all_residues)
                            
                            # Total atom count
                            atom_count = structure.numAtoms()
                            
                            # Additional properties from ProDy
                            if not properties:
                                properties = {}
                            properties['atom_count'] = atom_count
                            properties['chain_info'] = chain_info
                            if component_info:
                                properties['components'] = component_info
                            
                    finally:
                        # Clean up temp file if needed
                        if need_cleanup and os.path.exists(parse_path):
                            os.unlink(parse_path)
                    
            except Exception as e:
                # Fallback to simple parsing if ProDy fails
                if pdb_content:
                    chains = set()
                    residues = set()
                    for line in pdb_content.split('\n'):
                        if line.startswith(('ATOM', 'HETATM')):
                            if len(line) > 21:
                                chains.add(line[21])
                            if len(line) > 26:
                                try:
                                    residues.add((line[21], int(line[22:26].strip())))
                                except ValueError:
                                    pass
                    chain_count = len(chains)
                    residue_count = len(residues)
        
        protein = ProteinEntry(
            pdb_id=pdb_id,
            name=name,
            pdb_content=pdb_content,
            file_path=file_path,
            chain_count=chain_count,
            residue_count=residue_count,
            properties=properties or {},
            tags=tags or []
        )
        
        session = self.get_session(session_id)
        if session:
            session.proteins[protein_id] = protein
        
        return protein_id
    
    def get_protein(self, session_id: str, protein_id: str) -> Optional[ProteinEntry]:
        """Get protein from session"""
        session = self.get_session(session_id)
        return session.proteins.get(protein_id) if session else None
    
    def get_protein_by_pdb_id(self, session_id: str, pdb_id: str) -> Optional[tuple[str, ProteinEntry]]:
        """Get protein by PDB ID (returns protein_id and ProteinEntry)"""
        session = self.get_session(session_id)
        if not session:
            return None
        
        for protein_id, protein in session.proteins.items():
            if protein.pdb_id.lower() == pdb_id.lower():
                return (protein_id, protein)
        return None
    
    def list_proteins(self, session_id: str, tags: Optional[List[str]] = None) -> Dict[str, ProteinEntry]:
        """List proteins in session, optionally filtered by tags"""
        session = self.get_session(session_id)
        if not session:
            return {}
        
        if tags is None:
            return session.proteins.copy()
        
        filtered = {}
        for prot_id, protein in session.proteins.items():
            if any(tag in protein.tags for tag in tags):
                filtered[prot_id] = protein
        
        return filtered
    
    # File management methods
    def save_file(self, session_id: str, content: str, filename: str, file_type: str,
                 properties: Optional[Dict[str, Any]] = None, tags: Optional[List[str]] = None) -> str:
        """Save content to a file in the session's temp directory and return file ID"""
        session = self.get_session(session_id)
        if not session or not session.tmp_dir:
            raise ValueError(f"Session {session_id} not found or has no temp directory")
        
        file_id = str(uuid.uuid4())
        
        # Create safe filename with UUID to avoid collisions
        safe_filename = f"{file_id}_{filename}"
        file_path = os.path.join(session.tmp_dir, safe_filename)
        
        # Write content to file
        with open(file_path, 'w') as f:
            f.write(content)
        
        # Get file size
        size_bytes = os.path.getsize(file_path)
        
        # Create file entry
        file_entry = FileEntry(
            file_id=file_id,
            original_name=filename,
            file_path=file_path,
            file_type=file_type,
            size_bytes=size_bytes,
            properties=properties or {},
            tags=tags or []
        )
        
        session.files[file_id] = file_entry
        return file_id
    
    def get_file(self, session_id: str, file_id: str) -> Optional[FileEntry]:
        """Get file entry from session"""
        session = self.get_session(session_id)
        return session.files.get(file_id) if session else None
    
    def read_file_content(self, session_id: str, file_id: str) -> Optional[str]:
        """Read content of a file"""
        file_entry = self.get_file(session_id, file_id)
        if file_entry and os.path.exists(file_entry.file_path):
            with open(file_entry.file_path, 'r') as f:
                return f.read()
        return None
    
    def list_files(self, session_id: str, file_type: Optional[str] = None, 
                  tags: Optional[List[str]] = None) -> Dict[str, FileEntry]:
        """List files in session, optionally filtered by type or tags"""
        session = self.get_session(session_id)
        if not session:
            return {}
        
        files = session.files.copy()
        
        # Filter by file type if specified
        if file_type:
            files = {fid: f for fid, f in files.items() if f.file_type == file_type}
        
        # Filter by tags if specified
        if tags:
            files = {fid: f for fid, f in files.items() 
                    if any(tag in f.tags for tag in tags)}
        
        return files
    
    # Collection management
    def create_collection(self, session_id: str, collection_name: str, entity_ids: List[str]) -> bool:
        """Create a named collection of entities (molecules, proteins, or files)"""
        session = self.get_session(session_id)
        if session:
            session.collections[collection_name] = entity_ids
            return True
        return False
    
    def get_collection(self, session_id: str, collection_name: str) -> Dict[str, Any]:
        """Get entities from a named collection"""
        session = self.get_session(session_id)
        if not session or collection_name not in session.collections:
            return {}
        
        entity_ids = session.collections[collection_name]
        entities = {}
        
        for entity_id in entity_ids:
            # Check molecules
            if entity_id in session.molecules:
                entities[entity_id] = {'type': 'molecule', 'data': session.molecules[entity_id]}
            # Check proteins
            elif entity_id in session.proteins:
                entities[entity_id] = {'type': 'protein', 'data': session.proteins[entity_id]}
            # Check files
            elif entity_id in session.files:
                entities[entity_id] = {'type': 'file', 'data': session.files[entity_id]}
        
        return entities
    
    # Session management
    def clear_session(self, session_id: str) -> bool:
        """Clear all data for a session"""
        with self._lock:
            session = self._sessions.pop(session_id, None)
            if session and session.tmp_dir and os.path.exists(session.tmp_dir):
                try:
                    shutil.rmtree(session.tmp_dir)
                except Exception:
                    pass  # Best effort cleanup
            return session is not None
    
    def get_session_summary(self, session_id: str) -> Dict[str, Any]:
        """Get summary of session contents"""
        session = self.get_session(session_id)
        if not session:
            return {}
        
        # Calculate total file size
        total_file_size = sum(f.size_bytes for f in session.files.values())
        
        return {
            "session_id": session_id,
            "molecule_count": len(session.molecules),
            "protein_count": len(session.proteins),
            "file_count": len(session.files),
            "total_file_size_bytes": total_file_size,
            "total_file_size_mb": round(total_file_size / (1024 * 1024), 2),
            "collection_count": len(session.collections),
            "collections": list(session.collections.keys()),
            "tmp_directory": session.tmp_dir,
            "created_at": session.created_at.isoformat(),
            "last_accessed": session.last_accessed.isoformat()
        }

# Global session manager instance
session_manager = SessionStateManager()
