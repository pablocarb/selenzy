

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint, GetAtomPairFingerprint, GetTopologicalTorsionFingerprint
from rdkit.Chem.rdmolops import PatternFingerprint, RDKFingerprint
from rdkit import DataStructs
from os import path
import numpy as np

def getReactants(dbfile):
    clist = set()
    for line in open(dbfile):
        if line.startswith('#'):
            continue
        row= line.rstrip().split('\t')
        for x in row[1].split(' '):
            if x.startswith('MNXM'):
                clist.add(x)
    return clist

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

def getMols():
    mol = {}
    clist = getReactants(path.join('../data', 'reac_prop.tsv'))
    cstr = getStructs(path.join('../data', 'chem_prop.tsv'))
    for c in set(cstr) & clist:                  
        try:
            mol[c] = Chem.MolFromSmiles(cstr[c])
        except:
            continue
    return mol

def storeFingerprint(mol, ffun, fname, param=None, bit=False):
    fp = []
    fpNames = []
    print('Computing fingerprints...')
    for c in mol:
        try:
            if param is not None:
                fp.append( ffun(mol[c], param) )
            else:
                if bit:
                    fp.append( ffun(mol[c]).ToBitString()  )
                else:
                    fp.append( ffun(mol[c]) )
        except:
            continue
        fpNames.append(c)
    print('Saving...') 
    f = np.savez_compressed(fname, x=fp, y=fpNames)

def testPattern(ptfile, bit=False):
    """ Test how to reload PatternFingerprint """
    print('Validating fingerprint...')
    data = np.load(ptfile)
    fps = data['x']
    fpNames = data['y']
    if bit:
        fp = [DataStructs.CreateFromBitString(z) for z in fps]
    else:
        fp = fps
    sim = DataStructs.BulkTanimotoSimilarity(fp[0], list(fp))
    return fp


mol = getMols()
print('Pattern fingerprint....')
storeFingerprint(mol, PatternFingerprint, 'ptfp.npz', bit=True)
fp = testPattern('ptfp.npz', bit=True)          
print('RDK fingerprint....')
storeFingerprint(mol, RDKFingerprint, 'rdkfp.npz', bit=True)
fp = testPattern('rdkfp.npz', bit=True)
print('Atom pair fingerprint....')          
storeFingerprint(mol, GetAtomPairFingerprint, 'apfp.npz')
fp = testPattern('apfp.npz')
print('Topological torsion fingerprint....')          
storeFingerprint(mol, GetTopologicalTorsionFingerprint, 'ttfp.npz')
fp = testPattern('ttfp.npz')
print('Morgan fingerprint....')          
for radius in range(1,11):
     storeFingerprint(mol, GetMorganFingerprint, 'mgfp%d.npz' % (radius,) , radius)
     fp = testPattern('mgfp%d.npz' % (radius,))
