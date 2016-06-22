from rdkit import Chem
import rdkit.Chem.inchi
from urllib.error import HTTPError
import cirpy
import re

# In-Process conversions using rdkit:

def inchi_to_inchikey_get(inchi: str) -> str:
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    if inchikey == None:
        return ("Could not parse input: " + inchi, 500)

    return {"inchikey": inchikey}


def inchi_to_smiles_get(inchi: str) -> str:
    m = rdkit.Chem.inchi.MolFromInchi(inchi)
    if (m == None):
        return ("Could not parse input: " + inchi, 500)

    return {"smiles": Chem.MolToSmiles(m)}


def smiles_to_inchi_get(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    return {"inchi": rdkit.Chem.inchi.MolToInchi(m)}


def smiles_to_inchikey_get(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if (m == None):
        return ("Could not parse input: " + smiles, 500)

    inchi = rdkit.Chem.inchi.MolToInchi(m)
    inchikey = rdkit.Chem.inchi.InchiToInchiKey(inchi)
    return {"inchikey": inchikey}

# conversions using external REST services


class CirpyError(Exception):
    def __init__(self, code, message):
        self.code = code
        self.message = message


def resolve_via_cirpy(identifier: str, target: str, source: str,) -> str:
    try:
        converted = cirpy.resolve(identifier, target, [source])
        return converted
    except HTTPError as err:
        if err.code == 504 or err.code == 408:
            raise CirpyError(504, "Timeout while waiting for identifier resolution service")
        raise CirpyError(500, "HTTPError while communicating with identifier resolution service")


def inchikey_to_inchi_get(inchikey) -> str:
    try:
        inchi = resolve_via_cirpy(inchikey, 'stdinchi', 'stdinchikey')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)



def inchi_to_cas_get(inchi) -> str:
    try:
        cas = resolve_via_cirpy(inchi, 'cas', 'stdinchi')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_cas_get(inchikey) -> str:
    try:
        cas = resolve_via_cirpy(inchikey, 'cas', 'stdinchikey')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchi_get(cas) -> str:
    try:
        inchi = resolve_via_cirpy(cas, 'stdinchi', 'cas_number')
        return {"inchi": inchi}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_inchikey_get(cas) -> str:
    try:
        inchikey = clean_inchi_key(resolve_via_cirpy(cas, 'stdinchikey', 'cas_number'))
        return {"inchikey": inchikey}
    except CirpyError as err:
        return (err.message, err.code)


def cas_to_smiles_get(cas) -> str:
    try:
        smiles = resolve_via_cirpy(cas, 'smiles', 'cas_number')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)


def smiles_to_cas_get(smiles) -> str:
    try:
        cas = resolve_via_cirpy(smiles, 'cas', 'smiles')
        return {"cas": cas}
    except CirpyError as err:
        return (err.message, err.code)


def inchikey_to_smiles_get(inchikey) -> str:
    try:
        smiles = resolve_via_cirpy(inchikey, 'smiles', 'stdinchikey')
        return {"smiles": smiles}
    except CirpyError as err:
        return (err.message, err.code)

def clean_inchi_key(inchikey):
    return re.sub("(?i)InChiKey=", "", inchikey)
