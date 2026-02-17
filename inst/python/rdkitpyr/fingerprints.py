from rdkit import Chem, rdBase
from rdkit.Chem import MACCSkeys, rdmolops, rdFingerprintGenerator
from typing import Sequence
from .molecule_io import convert_to_molecules



def calculate_maccs_fps(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    verbose: bool = False,
):
    # Ref.: https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
    # Function: MACCSkeys.GenMACCSKeys()
    
    fp_size = 167 # hard-coded length for MACCS fingerprint
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        out = []
        molecules = convert_to_molecules(molecule_list)
        for molecule in molecules:
            if molecule is None:
                out.append([-1] * fp_size) # '-1' marks invalid molecules for R
            else:
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
    # Ref.: https://rdkit.org/docs/source/rdkit.Chem.rdmolops.html
    # Function: rdmolops.RDKFingerprint()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - atomInvariants
    #   - fromAtoms
    #   - atomBits
    #   - bitInfo
    
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        out = []
        molecules = convert_to_molecules(molecule_list)
        for molecule in molecules:
            if molecule is None:
                out.append([-1] * fpSize) # '-1' marks invalid molecules for R
            else:
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
    # Ref.: https://rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html
    # Function: generator.GetFingerprint()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - onlyNonzeroInvariants
    #   - countBounds
    #   - atomInvariantsGenerator
    #   - bondInvariantsGenerator
    #   - includeRedundantEnvironments
    
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
        
        out = []
        molecules = convert_to_molecules(molecule_list)
        for molecule in molecules:
            if molecule is None:
                out.append([-1] * fpSize) # '-1' marks invalid molecules for R
            else:
                fp = generator.GetFingerprint(molecule)
                out.append(list(fp))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out

