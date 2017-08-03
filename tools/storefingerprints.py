

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from os import path
import numpy as np

def getStructs(dbfile):
    structs = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        cid = m[0]
        smiles = m[6]
        if len(smiles) > 0:
            structs[cid] = smiles
    return structs

radius = 5
mol = {}
fp = []
fpNames = []
cstr = getStructs(path.join('data', 'chem_prop.tsv'))
for c in cstr:                  
    try:
        mol[c] = Chem.MolFromSmiles(cstr[c])
    except:
        continue
print('Computing fingerprints...')
for c in mol:
    try:
        fp.append( GetMorganFingerprint(mol[c], radius) )
    except:
        continue
    fpNames.append(c)
print('Saving...') 
f = np.savez('myfp3.npz', x=fp, y=fpNames)

