from rdkit import Chem
from rdkit.Chem import Descriptors



def calculate_all_descriptors(
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
                # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
                out.append(
                    Descriptors.CalcMolDescriptors(
                        molecule,
                        missingVal=None,
                        silent=not verbose
                    )
                )
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out



def calculate_exact_mass(
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
                # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
                out.append(Descriptors.ExactMolWt(molecule))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out



def calculate_molecular_weight(
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
                # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
                out.append(Descriptors.MolWt(molecule))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out


