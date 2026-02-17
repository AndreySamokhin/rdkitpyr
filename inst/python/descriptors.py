from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import Sequence



def calculate_all_descriptors(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    verbose: bool = False,
):
    # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
    # Descriptors.CalcMolDescriptors()
    #   Unused arguments:
    #     - silent
    
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    
    try:
        out = []
        empty_desc = None
        molecules = convert_to_molecules(molecule_list)
        for molecule in molecules:
            if molecule is None:
                # 'Descriptors.descList' is not part of public API
                if empty_desc is None:
                    dummy_mol = Chem.MolFromSmiles("C")
                    empty_desc = Descriptors.CalcMolDescriptors(dummy_mol)
                    empty_desc = {name: float('nan') for name in empty_desc}
                out.append(empty_desc.copy())
            else:
                out.append(
                    Descriptors.CalcMolDescriptors(
                        molecule,
                        missingVal=float('nan'),
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
        )
        for molecule in molecules:
            if molecule is None:
                out.append(float('nan'))
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
        )
        for molecule in molecules:
            if molecule is None:
                out.append(float('nan'))
            else:
                # https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
                out.append(Descriptors.MolWt(molecule))
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out


