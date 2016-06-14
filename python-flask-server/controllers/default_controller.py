from rdkit import Chem
import rdkit.Chem.inchi

# In-Process conversions using rdkit:

def inchi_to_inchikey_get(inchi) -> str:
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    if inchikey == None:
        return ("Could not parse input: " + inchi, 500)

    return {"inchikey": inchikey}


def inchi_to_smiles_get(inchi) -> str:
    m = rdkit.Chem.inchi.MolFromInchi(inchi)
    if (m == None):
        return ("Could not parse input: " + inchi, 500)

    return {"smiles": Chem.MolToSmiles(m)}


def smiles_to_inchi_get(smiles) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    return {"inchi": rdkit.Chem.inchi.MolToInchi(m)}


def smiles_to_inchikey_get(smiles) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    inchi = rdkit.Chem.inchi.MolToInchi(m)
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    return {"inchikey": inchikey}

# conversions using external REST services


def inchikey_to_inchi_get(inchikey) -> str:
    return 'do some magic!'


def inchi_to_cas_get(inchi) -> str:
    return 'do some magic!'


def inchikey_to_cas_get(inchikey) -> str:
    return 'do some magic!'


def cas_to_inchi_get(cas) -> str:
    return 'do some magic!'


def cas_to_inchikey_get(cas) -> str:
    return 'do some magic!'


def cas_to_smiles_get(cas) -> str:
    return 'do some magic!'


def smiles_to_cas_get(smiles) -> str:
    return 'do some magic!'


def inchikey_to_smiles_get(inchikey) -> str:
    return 'do some magic!'



