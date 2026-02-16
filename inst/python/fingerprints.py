from typing import Sequence
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from .molecule_io import convert_to_molecules


def generate_morgan_fingerprints(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    radius: int = 2,
    fp_size: int = 1024,
    count_simulation: bool = False,
    include_chirality: bool = False,
    use_bond_types: bool = True,
    only_nonzero_invariants: bool = False,
    include_ring_membership: bool = True,
    count_bounds: int | None = None,
):
    molecules = convert_to_molecules(molecule_list)
    ecfp_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        countSimulation=count_simulation,
        includeChirality=include_chirality,
        useBondTypes=use_bond_types,
        onlyNonzeroInvariants=only_nonzero_invariants,
        includeRingMembership=include_ring_membership,
        countBounds=count_bounds,
        fpSize=fp_size,
    )
    ecfp_fingerprints = []
    for mol in molecules:
        if mol is None:
            ecfp_fingerprints.append(np.full(fp_size, None))
        else:
            ecfp_fingerprints.append(ecfp_generator.GetFingerprintAsNumPy(mol))
    return np.vstack(ecfp_fingerprints)


def generate_rdkit_fingerprints(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    min_path: int = 1,
    max_path: int = 7,
    use_hs: bool = True,
    branched_paths: bool = True,
    use_bond_order: bool = True,
    count_simulation: bool = False,
    count_bounds: int | None = None,
    fp_size: int = 1024,
    num_bits_per_feature: int = 2,
):
    molecules = convert_to_molecules(molecule_list)
    rdkit_generator = rdFingerprintGenerator.GetRDKitFPGenerator(
        minPath=min_path,
        maxPath=max_path,
        useHs=use_hs,
        branchedPaths=branched_paths,
        useBondOrder=use_bond_order,
        countSimulation=count_simulation,
        countBounds=count_bounds,
        fpSize=fp_size,
        numBitsPerFeature=num_bits_per_feature,
    )
    rdkit_fingerprints = []
    for mol in molecules:
        if mol is None:
            rdkit_fingerprints.append(np.full(fp_size, None))
        else:
            rdkit_fingerprints.append(rdkit_generator.GetFingerprintAsNumPy(mol))
    return np.vstack(rdkit_fingerprints)


def generate_maccs_base(molecule):
    maccs = rdMolDescriptors.GetMACCSKeysFingerprint(molecule)
    maccs_arr = np.zeros(1)
    DataStructs.ConvertToNumpyArray(maccs, maccs_arr)
    return maccs_arr


def generate_maccs_keys(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
):
    molecules = convert_to_molecules(molecule_list)

    maccs_keys = []
    for mol in molecules:
        if mol is None:
            maccs_keys.append(np.full(167, None))
        else:
            maccs_keys.append(generate_maccs_base(mol))
    return np.vstack(maccs_keys)


def generate_atom_pairs_fingerprints(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    min_distance: int = 1,
    max_distance: int = 30,
    include_chirality: bool = False,
    use_2d: bool = True,
    count_simulation: bool = True,
    count_bounds: int | None = None,
    fp_size: int = 2048,
):
    molecules = convert_to_molecules(molecule_list)
    atom_pairs_generator = rdFingerprintGenerator.GetAtomPairGenerator(
        minDistance=min_distance,
        maxDistance=max_distance,
        includeChirality=include_chirality,
        use2D=use_2d,
        countSimulation=count_simulation,
        countBounds=count_bounds,
        fpSize=fp_size,
    )
    atom_pairs_fingerprints = []
    for mol in molecules:
        if mol is None:
            atom_pairs_fingerprints.append(np.full(fp_size, None))
        else:
            atom_pairs_fingerprints.append(
                atom_pairs_generator.GetFingerprintAsNumPy(mol)
            )
    return np.vstack(atom_pairs_fingerprints)


def generate_topological_torsion_fingerprints(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    torsion_atom_count: int = 4,
    include_chirality: bool = False,
    count_simulation: bool = True,
    count_bounds: int | None = None,
    fp_size: int = 2048,
):
    molecules = convert_to_molecules(molecule_list)
    topological_torsion_generator = (
        rdFingerprintGenerator.GetTopologicalTorsionGenerator(
            torsionAtomCount=torsion_atom_count,
            includeChirality=include_chirality,
            countSimulation=count_simulation,
            countBounds=count_bounds,
            fpSize=fp_size,
        )
    )
    topological_torsion_fingerprints = []
    for mol in molecules:
        if mol is None:
            topological_torsion_fingerprints.append(np.full(fp_size, None))
        else:
            topological_torsion_fingerprints.append(
                topological_torsion_generator.GetFingerprintAsNumPy(mol)
            )
    return np.vstack(topological_torsion_fingerprints)
