from rdkit import rdBase


def disable_rdkit_warnings():
    rdBase.DisableLog("rdApp.*") # disable RDKit warning messages globally


def enable_rdkit_warnings():
    rdBase.EnableLog("rdApp.*") # enable RDKit warning messages globally
    
    