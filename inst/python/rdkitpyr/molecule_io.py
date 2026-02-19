from rdkit import Chem
from typing import Sequence, Optional


def parse_molecules(
    molecule_strings: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
):
    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.inchi.html
    # Function: inchi.MolFromInchi()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - sanitize
    #   - removeHs
    #   - logLevel
    #   - treatWarningAsError

    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html
    # Function: rdmolfiles.MolFromSmiles()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - sanitize
    #   - replacements
    
    if isinstance(molecule_strings, str) or isinstance(molecule_strings, Chem.Mol):
        molecule_strings = [molecule_strings]
    
    if len(molecule_strings) == 0:
        return []
        
    out = []
    for molecule in molecule_strings:
        try:
            if isinstance(molecule, str):
                if molecule.startswith("InChI="):
                    mol = Chem.MolFromInchi(molecule)
                else:
                    mol = Chem.MolFromSmiles(molecule)
                out.append(mol)
            elif isinstance(molecule, Chem.Mol):
                out.append(molecule)
            else:
                out.append(None)
        except Exception:
            out.append(None)
    return out


def convert_to_smiles(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
):
    # Ref.: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html
    # Function: rdmolfiles.MolToSmiles()
    # Note: The following arguments exist but are not used in this wrapper:
    #   - rootedAtAtom
    #   - doRandom
    #   - ignoreAtomMapNumbers 
    
    out = []
    molecules = parse_molecules(molecule_list)
    for molecule in molecules:
        if molecule is None:
            out.append("")
        else:
            try:
                out.append(
                    Chem.MolToSmiles(
                        molecule,
                        isomericSmiles=isomericSmiles,
                        kekuleSmiles=kekuleSmiles,
                        canonical=canonical,
                        allBondsExplicit=allBondsExplicit,
                        allHsExplicit=allHsExplicit,
                    )
                )
            except Exception:
                out.append("")
    return out


def convert_to_inchi(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    removeHs: bool = True,
    replacements: Optional[dict] = None,
):
    molecules = parse_molecules(
        molecule_list,
        sanitize=sanitize,
        removeHs=removeHs,
        replacements=replacements,
    )
    
    out = []
    for molecule in molecules:
        if molecule is None:
            out.append("")
        else:
            try:
                out.append(Chem.MolToInchi(molecule))
            except Exception:
                out.append("")
    return out


def convert_to_inchikey(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    replacements: Optional[dict] = None,
):
    if replacements is None:
        replacements = {}
    if isinstance(molecule_list, str) or isinstance(molecule_list, Chem.Mol):
        molecule_list = [molecule_list]
        
    if len(molecule_list) == 0:
        return []
    
    out = []
    for molecule in molecule_list:
        try:
            if isinstance(molecule, str):
                if molecule.startswith("InChI="):
                    inchikey = Chem.InchiToInchiKey(molecule)
                    if inchikey is None:
                        out.append("")
                    else:
                        out.append(inchikey)
                else:
                    mol = Chem.MolFromSmiles(
                        molecule, sanitize=sanitize, replacements=replacements
                    )
                    if mol is not None:
                        out.append(Chem.MolToInchiKey(mol))
                    else:
                        out.append("")
            elif isinstance(molecule, Chem.Mol):
                out.append(Chem.MolToInchiKey(molecule))
            else:
                out.append("")
        except Exception:
            out.append("")
    return out

