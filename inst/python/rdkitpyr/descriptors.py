from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors
from typing import Sequence
from .molecule_io import parse_molecules



def calculate_all_descriptors(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
):
    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
    # Function: Descriptors.CalcMolDescriptors()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - silent
    
    out = []
    empty_desc = None
    molecules = parse_molecules(molecule_list)
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
    return out


def calculate_exact_mass(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
):
    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
    # Function: Descriptors.ExactMolWt()
    
    out = []
    molecules = parse_molecules(molecule_list)
    for molecule in molecules:
        if molecule is None:
            out.append(float('nan'))
        else:
            out.append(Descriptors.ExactMolWt(molecule))
    return out



def calculate_molecular_weight(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
):
    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
    # Function: Descriptors.MolWt()
    
    out = []
    molecules = parse_molecules(molecule_list)
    for molecule in molecules:
        if molecule is None:
            out.append(float('nan'))
        else:
            out.append(Descriptors.MolWt(molecule))
    return out


