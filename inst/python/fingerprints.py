from rdkit import Chem
from rdkit.Chem import MACCSkeys, rdmolops, rdFingerprintGenerator



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
            verbose=verbose,
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


