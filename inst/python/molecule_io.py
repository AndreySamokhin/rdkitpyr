from rdkit import Chem
from rdkit import rdBase
from typing import Sequence, Optional


def convert_to_molecules(
    molecule_strings: str | Chem.Mol | Sequence[str] | Sequence[Chem.Mol],
    sanitize: bool = True,
    removeHs: bool = True,
    replacements: Optional[dict] = None,
    verbose: bool = False,
):
    out = []
    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:
        if isinstance(molecule_strings, str):
            molecule_strings = [molecule_strings]
        if len(molecule_strings) == 0:
            return []

        for molecule in molecule_strings:
            if isinstance(molecule, str):
                if molecule.startswith("InChI="):
                    mol = Chem.MolFromInchi(
                        molecule, sanitize=sanitize, removeHs=removeHs
                    )
                else:
                    if replacements is None:
                        replacements = {}
                    mol = Chem.MolFromSmiles(
                        molecule, sanitize=sanitize, replacements=replacements
                    )
                out.append(mol)
            elif isinstance(molecule, Chem.Mol):
                out.append(molecule)
            else:
                out.append(None)
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
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
            sanitize=sanitize,
            removeHs=removeHs,
            replacements=replacements,
            verbose=verbose,
        )
        for molecule in molecules:
            if molecule is None:
                out.append(None)
            else:
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
            sanitize=sanitize,
            removeHs=removeHs,
            replacements=replacements,
            verbose=verbose,
        )
        for molecule in molecules:
            if molecule is None:
                out.append(None)
            else:
                out.append(
                    Chem.MolToInchi(
                        molecule,
                    )
                )
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
    out = []

    if not verbose:
        rdBase.DisableLog("rdApp.*")
    try:

        if isinstance(molecule_list, str) or isinstance(molecule_list, Chem.Mol):
            molecule_list = [molecule_list]
        if len(molecule_list) == 0:
            return []
        for molecule in molecule_list:
            if isinstance(molecule, str):
                if molecule.startswith("InChI="):
                    out.append(Chem.InchiToInchiKey(molecule))
                else:
                    if replacements is None:
                        replacements = {}
                    mol = Chem.MolFromSmiles(
                        molecule, sanitize=sanitize, replacements=replacements
                    )
                    if mol is not None:
                        inchikey = Chem.MolToInchiKey(mol)
                        out.append(inchikey)
                    else:
                        out.append(None)
            elif isinstance(molecule, Chem.Mol):
                out.append(Chem.MolToInchiKey(molecule))
            else:
                out.append(None)
    finally:
        if not verbose:
            rdBase.EnableLog("rdApp.*")
    return out
