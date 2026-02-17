from rdkit import Chem, rdBase
from typing import Sequence, Optional


def convert_to_molecules(
    molecule_strings: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    removeHs: bool = True,
    replacements: Optional[dict] = None
):
    if replacements is None:
        replacements = {}

    if isinstance(molecule_strings, str) or isinstance(molecule_strings, Chem.Mol):
        molecule_strings = [molecule_strings]
    
    if len(molecule_strings) == 0:
        return []
        
    out = []
    for molecule in molecule_strings:
        try:
            if isinstance(molecule, str):
                if molecule.startswith("InChI="):
                    mol = Chem.MolFromInchi(
                        molecule, sanitize=sanitize, removeHs=removeHs
                    )
                else:
                    mol = Chem.MolFromSmiles(
                        molecule, sanitize=sanitize, replacements=replacements
                    )
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
    sanitize: bool = True,
    removeHs: bool = True,
    replacements: Optional[dict] = None,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    rootedAtAtom: int = -1,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
    doRandom: bool = False,
    ignoreAtomMapNumbers: bool = False,
    verbose: bool = False,
):
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        molecules = convert_to_molecules(
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
                    out.append(
                        Chem.MolToSmiles(
                            molecule,
                            isomericSmiles=isomericSmiles,
                            kekuleSmiles=kekuleSmiles,
                            rootedAtAtom=rootedAtAtom,
                            canonical=canonical,
                            allBondsExplicit=allBondsExplicit,
                            allHsExplicit=allHsExplicit,
                            doRandom=doRandom,
                            ignoreAtomMapNumbers=ignoreAtomMapNumbers,
                        )
                    )
                except Exception:
                    out.append("")
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out


def convert_to_inchi(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    removeHs: bool = True,
    replacements: Optional[dict] = None,
    verbose: bool = False,
):
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        molecules = convert_to_molecules(
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
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out


def convert_to_inchikey(
    molecule_list: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    replacements: Optional[dict] = None,
    verbose: bool = False,
):
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    
    if replacements is None:
        replacements = {}
    
    try:
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
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out

