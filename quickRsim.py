'''
quickRsim (c) University of Manchester 2017

quickRsim is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  Pablo Carbonell, SYNBIOCHEM
@description: Compute a fast reaction similarity with the database
@examples: 
1. Compute similarity to one reaction in the database (requires -chem parameter)
python quickRsim.py data/reac_prop.tsv data/metanetx.fs -rid MNXR3215 -chem data/chem_prop.tsv
2. Compute similarity to a given reaction file for a threshold above 0.9
python quickRsim.py data/reac_prop.tsv data/metanetx.fs -rxn rhea15870.rxn -th 0.9
'''

from __future__ import print_function
import argparse
import subprocess
import re
from os import path
import os
from rdkit.Chem import rdChemReactions
import math
import numpy as np
from rdkit import DataStructs, Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint, GetConnectivityInvariants
from rdkit.Chem import AllChem


def storeReaction(smi, rfile):
    left, right = smi.split('>>')
    subs = left.split('.')
    prods = right.split('.')
    sd, pd = ({},{})
    for s in subs:
        if s not in sd:
            sd[s] = 0
        sd[s] += 1
    for p in prods:
        if p not in pd:
            pd[p] = 0
        pd[p] += 1
    rsp = {rfile: (sd, pd)}
    return rsp

def getReaction(rfile):
    rxn = rdChemReactions.ReactionFromRxnFile(rfile)
    smi = rdChemReactions.ReactionToSmiles(rxn)
    return storeReaction(smi, rfile)

def getReactionFromSmiles(smi, rxnfile):
    rxn = rdChemReactions.ReactionFromSmarts(smi)
    mdl = rdChemReactions.ReactionToRxnBlock(rxn)
    with open(rxnfile, 'w') as handler:
        handler.write(mdl)
    return storeReaction(smi, rxnfile)

def getReactionFromSmilesFile(smartsfile, rxnfile):
    with open(smartsfile) as handler:
        smarts = handler.readline()
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    smi = rdChemReactions.ReactionToSmiles(rxn)
    mdl = rdChemReactions.ReactionToRxnBlock(rxn)
    with open(rxnfile, 'w') as handler:
        handler.write(mdl)
    return storeReaction(smi, rxnfile)

def getClosest(smi, fpfile, th=0.8):
    dist = {}
    data = np.load(fpfile)
    fp = data['x'] 
    fpNames = data['y']
    data.close()
    radius = 5
    targetMol = Chem.MolFromSmiles(smi)
    # If RDkit fails, we sanitize first using molconvert from ChemAxon, which is more robust
    if targetMol is None:
        cmd = ['molconvert', 'mol', smi]
        cmd2 = ['molconvert', 'smiles']
        job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        job2 = subprocess.Popen(cmd2, stdin=job.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        job.stdout.close()
        out, err = job2.communicate()
        targetMol = Chem.MolFromSmiles(out)

    targetFp = GetMorganFingerprint(targetMol, radius)
    tn = DataStructs.BulkTanimotoSimilarity(targetFp, list(fp))
    for i in sorted(range(0, len(tn))):
        dist[fpNames[i]] = tn[i]
    return dist


def getReactants(equation):
    reactants = {}
    for x in equation.split(' + '):
        n, c = x.split(' ')
        try:
            n = int(n)
        except:
            pass
        reactants[c] = n
    return reactants

def reacSubsProds(dbfile):
    rsp = {}
    for l in open(dbfile):
        if l.startswith('#'):
            continue
        m = l.rstrip().split('\t')
        rid = m[0]
        subs = {}
        prods = {}
        m = l.rstrip().split('\t')
        left, right = m[1].split(' = ')
        subs = getReactants(left)
        prods = getReactants(right)
        if len(subs) > 0 and len(prods) > 0:
            rsp[rid] = (subs, prods)
    return rsp

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

def getRSim(s1, p1, s2, p2, sim):
    cl = {'s1': s1, 'p1': p1, 's2': s2, 'p2':p2}
    ss = {} 
    simm = {}
    pairs = [('s1','s2'), ('s1', 'p2'), ('p1', 's2'), ('p1', 'p2')]
    for p in pairs:
        pairings = set()
        simm[p] = {}
        for x in cl[p[0]]:
            simm[p][x] = (0.0, x, None)
            if x in sim:
                for y in cl[p[1]]:
                    if y in sim[x]:
                        pairings.add( (sim[x][y], x, y) )
                           # import pdb
                            #pdb.set_trace()
        found = {'left': set(), 'right': set()}
        for v in sorted(pairings, key = lambda h: -h[0]):
            if v[1] not in found['left'] and v[2] not in found['right']:
                if v[0] > simm[p][v[1]][0]:
                    simm[p][v[1]] = v
                    found['left'].add(v[1])
                    found['right'].add(v[2])
        s = []
        for x in simm[p]:
            s.append(simm[p][x][0])
        if len(s) > 0:
            ss[p] = sum(s)/len(s)
        else:
            ss[p] = 0.0
    S1 = math.sqrt(ss[pairs[0]]**2 + ss[pairs[3]]**2)/math.sqrt(2)
    S2 = math.sqrt(ss[pairs[1]]**2 + ss[pairs[2]]**2)/math.sqrt(2)
    return(S1, S2)

def arguments():
    parser = argparse.ArgumentParser(description='quickRSim Pablo Carbonell, SYNBIOCHEM, 2016')
    parser.add_argument('db', help='Metanetx reaction file')
    parser.add_argument('fp', help='Babel fingerprint file for reactants')
    parser.add_argument('-rxn', 
                        help='Input reaction rxn file')
    parser.add_argument('-rid', 
                        help='Input reaction id')
    parser.add_argument('-smarts', 
                        help='Input reaction SMARTS')
    parser.add_argument('-smartsfile', 
                        help='Input reaction SMARTS file')
    parser.add_argument('-chem', 
                        help='Metanetx chemical structures (if input is reaction id)')
    parser.add_argument('-th', type=float, default=0.8, 
                        help='Similarity threshold [default=0.8]')
    parser.add_argument('-out', 
                        help='Output results in .txt file, please specify file name')
    parser.add_argument('-high', 
                        help='Output results in .txt file with highest similarity score from both forwards and backwards reactions, please specify file name')
    arg = parser.parse_args()
    return arg



if __name__ == '__main__':
    arg = arguments()
    if arg.out:
        fileObj = open(arg.out, 'w')
    
    if arg.high:
        fileObj = open(arg.high, 'w')
    rsp = reacSubsProds(arg.db)
    if arg.rxn is not None:
        rTarget = getReaction(arg.rxn)
    elif arg.smarts is not None:
        rxnfile = os.path.join(os.path.dirname(arg.out), 'reaction.rxn')
        rTarget = getReactionFromSmiles(arg.smarts, rxnfile)
    elif arg.smartsfile is not None:
        rxnfile = os.path.join(os.path.dirname(arg.out), 'reaction.rxn')
        rTarget = getReactionFromSmilesFile(arg.smartsfile, rxnfile)        
    elif arg.rid is not None:
        struct = getStructs(arg.chem)
        rTarget = {arg.rid: [{},{}]}
        for side in (0,1):
            for s in rsp[arg.rid][side]:
                if s in struct:
                    rTarget[arg.rid][side][struct[s]] = rsp[arg.rid][side][s]
    else:
        raise Exception('No target')
    sim = {} 
    # Compute similarities for new reactants
    for r in rTarget:
        for s in rTarget[r][0]:
            if s not in sim:
                try:
                    sim[s] = getClosest(s, arg.fp, arg.th)
                except:
                    continue
        for p in rTarget[r][1]:
            if p not in sim:
                try:
                    sim[p] = getClosest(p, arg.fp, arg.th)
                except:
                    continue
    # Get reaction similarities
    for r1 in rTarget:
        s1, p1 = rTarget[r1]
        for r2 in rsp:
#            if r2 == 'MNXR70380':
#                import pdb
#                pdb.set_trace()
            s2, p2 = rsp[r2]
            debug = False
            S1, S2 = getRSim(s1, p1, s2, p2, sim)
            if S1 > 0 and S2 > 0:
                print(r1, r2, S1, S2)
                
                if arg.out:
                    print(r1, r2, S1, S2, file = fileObj)
                
                if arg.high:
                    if S1 >= S2:
                        print(r1, r2, S1, file = fileObj)
                    else:
                        print(r1, r2, S2, file = fileObj)
                        

                        
