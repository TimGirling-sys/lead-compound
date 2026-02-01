"""
Fingerprint calculation module for molecular structures.

This module handles SDF file parsing and molecular fingerprint generation
using RDKit's Morgan/ECFP fingerprints.
"""

from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass, field
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit import DataStructs


@dataclass
class Compound:
    """
    Represents a chemical compound with its molecular data.

    Attributes:
        mol: RDKit molecule object
        name: Compound name/identifier
        smiles: SMILES representation
        fingerprint: Morgan fingerprint bit vector
        properties: Dictionary of additional properties from SDF
        mol_weight: Molecular weight
        num_heavy_atoms: Number of heavy (non-hydrogen) atoms
        index: Position in the original SDF file
    """
    mol: Chem.Mol
    name: str
    smiles: str
    fingerprint: DataStructs.ExplicitBitVect
    properties: Dict[str, Any] = field(default_factory=dict)
    mol_weight: float = 0.0
    num_heavy_atoms: int = 0
    index: int = 0

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other):
        if isinstance(other, Compound):
            return self.smiles == other.smiles
        return False


class FingerprintCalculator:
    """
    Calculates molecular fingerprints and parses SDF files.

    Uses Morgan (ECFP) fingerprints for structural similarity calculations.
    Morgan fingerprints are circular fingerprints that capture local
    structural features around each atom.

    Attributes:
        radius: Radius for Morgan fingerprint calculation (default 2 = ECFP4)
        n_bits: Number of bits in fingerprint vector (default 2048)
    """

    def __init__(self, radius: int = 2, n_bits: int = 2048):
        """
        Initialize the fingerprint calculator.

        Args:
            radius: Morgan fingerprint radius.
                   2 = ECFP4 (captures up to 4-bond environments)
                   3 = ECFP6 (captures up to 6-bond environments)
            n_bits: Length of the bit vector for fingerprints.
                   Higher values reduce collision probability.
        """
        self.radius = radius
        self.n_bits = n_bits

    def parse_sdf(self, sdf_path: str) -> List[Compound]:
        """
        Parse an SDF file and extract all valid compounds.

        Args:
            sdf_path: Path to the SDF file

        Returns:
            List of Compound objects with fingerprints calculated

        Raises:
            FileNotFoundError: If SDF file doesn't exist
            ValueError: If no valid compounds found
        """
        supplier = Chem.SDMolSupplier(sdf_path, removeHs=True)
        compounds = []

        for idx, mol in enumerate(supplier):
            if mol is None:
                continue

            compound = self._mol_to_compound(mol, idx)
            if compound is not None:
                compounds.append(compound)

        if not compounds:
            raise ValueError(f"No valid compounds found in {sdf_path}")

        return compounds

    def parse_smiles_list(self, smiles_list: List[str],
                          names: Optional[List[str]] = None) -> List[Compound]:
        """
        Parse a list of SMILES strings into Compound objects.

        Args:
            smiles_list: List of SMILES strings
            names: Optional list of compound names

        Returns:
            List of Compound objects with fingerprints calculated
        """
        if names is None:
            names = [f"Compound_{i}" for i in range(len(smiles_list))]

        compounds = []
        for idx, (smi, name) in enumerate(zip(smiles_list, names)):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                compound = self._mol_to_compound(mol, idx)
                if compound is not None:
                    compound.name = name
                    compounds.append(compound)

        return compounds

    def _mol_to_compound(self, mol: Chem.Mol, index: int) -> Optional[Compound]:
        """
        Convert an RDKit molecule to a Compound object.

        Args:
            mol: RDKit molecule object
            index: Index in the source file/list

        Returns:
            Compound object or None if conversion fails
        """
        try:
            # Calculate fingerprint
            fp = GetMorganFingerprintAsBitVect(
                mol,
                radius=self.radius,
                nBits=self.n_bits
            )

            # Get SMILES
            smiles = Chem.MolToSmiles(mol, canonical=True)

            # Get name from molecule properties
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Mol_{index}"
            if not name or name.strip() == "":
                name = f"Mol_{index}"

            # Extract all properties
            properties = {}
            for prop_name in mol.GetPropsAsDict():
                if not prop_name.startswith("_"):
                    properties[prop_name] = mol.GetProp(prop_name)

            # Calculate molecular properties
            mol_weight = Descriptors.MolWt(mol)
            num_heavy_atoms = mol.GetNumHeavyAtoms()

            return Compound(
                mol=mol,
                name=name,
                smiles=smiles,
                fingerprint=fp,
                properties=properties,
                mol_weight=mol_weight,
                num_heavy_atoms=num_heavy_atoms,
                index=index
            )

        except Exception as e:
            print(f"Warning: Failed to process molecule at index {index}: {e}")
            return None

    def calculate_similarity(self, compound1: Compound,
                            compound2: Compound) -> float:
        """
        Calculate Tanimoto similarity between two compounds.

        Tanimoto coefficient = |A ∩ B| / |A ∪ B|
        where A and B are the fingerprint bit sets.

        Args:
            compound1: First compound
            compound2: Second compound

        Returns:
            Tanimoto similarity coefficient (0.0 to 1.0)
        """
        return DataStructs.TanimotoSimilarity(
            compound1.fingerprint,
            compound2.fingerprint
        )

    def calculate_similarity_matrix(self,
                                    compounds: List[Compound]) -> np.ndarray:
        """
        Calculate pairwise Tanimoto similarity matrix for all compounds.

        Args:
            compounds: List of compounds

        Returns:
            NxN numpy array of similarity values
        """
        n = len(compounds)
        similarity_matrix = np.zeros((n, n))

        for i in range(n):
            similarity_matrix[i, i] = 1.0
            for j in range(i + 1, n):
                sim = self.calculate_similarity(compounds[i], compounds[j])
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim

        return similarity_matrix

    def get_similar_compounds(self,
                             target: Compound,
                             compounds: List[Compound],
                             threshold: float = 0.7) -> List[Tuple[Compound, float]]:
        """
        Find compounds similar to a target compound.

        Args:
            target: Target compound to compare against
            compounds: List of compounds to search
            threshold: Minimum similarity threshold

        Returns:
            List of (compound, similarity) tuples sorted by similarity
        """
        similar = []
        for compound in compounds:
            if compound.smiles != target.smiles:
                sim = self.calculate_similarity(target, compound)
                if sim >= threshold:
                    similar.append((compound, sim))

        return sorted(similar, key=lambda x: x[1], reverse=True)
