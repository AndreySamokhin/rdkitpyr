from rdkit import Chem
from rdkit.Chem import MACCSkeys, rdmolops, rdFingerprintGenerator
from typing import Sequence



def calculate_maccs_fps(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    verbose: bool = False,
):
    out = []

    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        if isinstance(molecule_list, str) or isinstance(molecule_list, Chem.Mol):
            molecule_list = [molecule_list]
        if len(molecule_list) == 0:
            return []
        molecules = convert_to_molecules(
            molecule_list,
        )
        for molecule in molecules:
            if molecule is None:
                out.append(None)
            else:
                # https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
                fp = MACCSkeys.GenMACCSKeys(molecule)
                out.append(list(fp))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out



def calculate_rdkit_fps(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    minPath: int = 1,
    maxPath: int = 7,
    fpSize: int = 2048,
    nBitsPerHash: int = 2,
    useHs: bool = True,
    tgtDensity: float = 0.0,
    minSize: int = 128,
    branchedPaths: bool = True,
    useBondOrder: bool = True,
    verbose: bool = False,
):
    # TODO: The 'RDKFingerprint()' function contains some additional arguments:
    #   - atomInvariants
    #   - fromAtoms
    #   - atomBits
    #   - bitInfo
    
    out = []

    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        if isinstance(molecule_list, str) or isinstance(molecule_list, Chem.Mol):
            molecule_list = [molecule_list]
        if len(molecule_list) == 0:
            return []
        molecules = convert_to_molecules(
            molecule_list,
        )
        for molecule in molecules:
            if molecule is None:
                out.append(None)
            else:
                # https://rdkit.org/docs/source/rdkit.Chem.rdmolops.html
                fp = rdmolops.RDKFingerprint(
                    molecule,
                    minPath=minPath,
                    maxPath=maxPath,
                    fpSize=fpSize,
                    nBitsPerHash=nBitsPerHash,
                    useHs=useHs,
                    tgtDensity=tgtDensity,
                    minSize=minSize,
                    branchedPaths=branchedPaths,
                    useBondOrder=useBondOrder
                )
                out.append(list(fp))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out



def calculate_morgan_fps(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    radius: int = 3,
    fpSize: int = 2048,
    countSimulation: bool = False,
    includeChirality: bool = False,
    useBondTypes: bool = True,
    includeRingMembership: bool = True,
    verbose: bool = False,
):
    # TODO: The 'GetMorganGenerator()' function contains some additional
    # arguments:
    #   - onlyNonzeroInvariants
    #   - countBounds
    #   - atomInvariantsGenerator
    #   - bondInvariantsGenerator
    #   - includeRedundantEnvironments
  
    out = []

    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius,
            fpSize=fpSize,
            countSimulation=countSimulation,
            includeChirality=includeChirality,
            useBondTypes=useBondTypes,
            includeRingMembership=includeRingMembership
        )
        if isinstance(molecule_list, str) or isinstance(molecule_list, Chem.Mol):
            molecule_list = [molecule_list]
        if len(molecule_list) == 0:
            return []
        molecules = convert_to_molecules(
            molecule_list,
        )
        for molecule in molecules:
            if molecule is None:
                out.append(None)
            else:
                # https://rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html
                fp = generator.GetFingerprint(molecule)
                out.append(list(fp))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out


