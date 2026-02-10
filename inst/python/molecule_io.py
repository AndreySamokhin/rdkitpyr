from rdkit import Chem


def convert_inchi_to_smiles(inchi_list, **kwargs):
    # **kwargs - Additional keyword arguments passed to Chem.MolToSmiles()
    #            (e.g., isomericSmiles=True, kekuleSmiles=False, canonical=True)
    out = []
    for inchi in inchi_list:
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                out.append(None)
            else:
                out.append(Chem.MolToSmiles(mol, **kwargs))
        except:
            out.append(None)
    return out


def convert_inchi_to_inchikey(inchi_list):
    out = []
    for inchi in inchi_list:
        try:
            out.append(Chem.InchiToInchiKey(inchi))
        except:
            out.append(None)
    return out


def convert_inchi_to_inchi(inchi_list):
    out = []
    for inchi in inchi_list:
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                out.append(None)
            else:
                out.append(Chem.MolToInchi(mol))
        except:
            out.append(None)
    return out


def convert_smiles_to_inchi(smiles_list):
    out = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out.append(None)
            else:
                out.append(Chem.MolToInchi(mol))
        except:
            out.append(None)
    return out


def convert_smiles_to_inchikey(smiles_list):
    out = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                out.append(None)
            else:
                inchi = Chem.MolToInchi(mol)
                out.append(Chem.InchiToInchiKey(inchi))
        except:
            out.append(None)
    return out


def convert_smiles_to_smiles(smiles_list, **kwargs):
    # **kwargs - Additional keyword arguments passed to Chem.MolToSmiles()
    #            (e.g., isomericSmiles=True, kekuleSmiles=False, canonical=True)
    out = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                out.append(None)
            else:
                out.append(Chem.MolToSmiles(mol, **kwargs))
        except:
            out.append(None)
    return out
    
