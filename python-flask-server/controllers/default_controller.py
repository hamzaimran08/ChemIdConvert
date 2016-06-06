from rdkit import Chem
from rdkit.Chem import inchi

def smiles_to_inchi_get(smiles) -> str:
    m = Chem.MolFromSmiles(smiles)
    return inchi.MolToInchi(m)
