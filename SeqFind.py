#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 10:53:03 2017

@author: jerrywzy
"""
import re
import os, subprocess
import json
import csv
import argparse

def readFasta(fileFasta):
    
    from Bio import SeqIO
    
    sequence = {}
    descriptions = {}
    osource = {}
    names = []
    seen = set()
    seen_add = seen.add
    
    for seq_record in SeqIO.parse(fileFasta, "fasta"):
        ming = seq_record.id
        idonly = re.search(r'\|(.*?)\|',ming)
        x = idonly.group(1)
        
        fulldesc = seq_record.description
        desc = fulldesc.rsplit('OS=')[0]
        shortdesc = " ".join(desc.split()[1:])
        orgsource = fulldesc.rsplit('OS=')[1]
        shortos = orgsource.rsplit('GN=')[0]
    #    o = shortos.group(1)
        if ',' in shortdesc:
            y = shortdesc.replace(",", ";")
        else:
            y = shortdesc
        
        if x not in seen:
            names.append(x)
            seen_add(x)
        myseq = seq_record.seq
        sequence[x]=(myseq)
        descriptions[x]= y
        osource[x]=shortos
    
        
#    with open('sequences.json', 'w') as f:
#        json.dump(sequence, f)
    
#    with open('descriptions.json', 'w') as f2:
#        json.dump(descriptions, f2)
    
    
    return (sequence, names, descriptions, osource)

def readReacFile(reacfile):
    fileObj2 = open(reacfile, 'r')
    MnxToUprot = dict()
    for line in fileObj2:
        splitdata = line.split()
        Mnx = splitdata[0]
        Uprot = splitdata[2]
        MnxToUprot.setdefault(Mnx, []).append(Uprot)
    fileObj2.close()
    
    with open('MnxToUprot.json', 'w') as f:
        json.dump(MnxToUprot, f)
    
    return (MnxToUprot)
 
def readRxnCons(consensus):
    
    f = open(consensus, 'r')
    
    MnxDir = {}
    
    for line in f:
        splitdata = line.split()
        Mnx = splitdata[1]
        dirxn = splitdata[2]
        MnxDir[Mnx] = dirxn
        
    return (MnxDir)   
    
def getMnxSim(rxn, drxn=0):
    cmd = ['/home/jerrywzy/anaconda3/envs/p2/bin/python', 'quickRsim.py', 
           'data/reac_prop.tsv', 'data/fp.npz', '-rxn', rxn, '-out', 'results_rsim.txt']
    job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = job.communicate()
#    print(out)
#    global mnx
#    mnx = []
    MnxSim = {}
    MnxDirPref = readRxnCons("rxn_consensus_20160612.txt")
    MnxDirUsed = {}
    
    if drxn==1:
        file = open("results_rsim.txt", 'r')
        for line in file:
            splitdata = line.split()
            S1 = splitdata[2]
            S2 = splitdata[3]
            Mnx = splitdata[1]

            try:
                direction = MnxDirPref[Mnx] 
                if direction == '1':
                    MnxSim[Mnx]=S1
                    MnxDirUsed[Mnx]= '1'
                elif direction == '-1':
                    MnxSim[Mnx]=S2
                    MnxDirUsed[Mnx]= '-1'
                else:
                    if S1 >= S2:
                        MnxSim[Mnx]=S1
                        MnxDirUsed[Mnx]= '1'
                    elif S2 > S1:
                        MnxSim[Mnx]=S2
                        MnxDirUsed[Mnx]= '-1'
            except KeyError:
                MnxDirPref[Mnx] = "N/A"   # for missing Keys?
                if S1 >= S2:
                    MnxSim[Mnx]=S1
                else:
                    MnxSim[Mnx]=S2
        file.close()
            
        return (MnxSim, MnxDirPref, MnxDirUsed)
        
    else:
        
        file = open("results_rsim.txt", 'r')
        for line in file:
            splitdata = line.split()
            S1 = splitdata[2]
            S2 = splitdata[3]
            Mnx = splitdata[1]
            if S1 > S2:
                MnxDirUsed[Mnx]='1'
                MnxSim[Mnx]=S1
            elif S2 > S1:
                MnxDirUsed[Mnx]='-1'
                MnxSim[Mnx]=S2
            elif S2 == S1:
                MnxDirUsed[Mnx]='0'
                MnxSim[Mnx]=S2                
        file.close()              
        
        return (MnxSim, MnxDirPref, MnxDirUsed)

   
def pepstats(file):
    args = ("pepstats -sequence {0} -outfile results.pepstats".format(file))
    os.system(args)
    
    f = open("results.pepstats", "r")
    
    hydrop = {}
    weight = {}
    isoelec = {}
    polar = {}
    
    for line in f:
        if "PEPSTATS of" in line:
            splitdata = line.split()
            seq = splitdata[2]
        elif "Molecular weight = " in line:
            splitdata = line.split()
            w = splitdata[3]
            weight[seq] = w
        elif "Isoelectric Point = " in line:
            splitdata = line.split()
            i = splitdata[3]
            isoelec[seq] = i
        elif "Polar	" in line:
            splitdata = line.split()
            percent = splitdata[3]
            polar[seq] = percent
            
    return (hydrop, weight, isoelec, polar)
        
def garnier(file):
    args = ("garnier -sequence {0} -outfile garnier.txt".format(file))
    os.system(args)
    
    f = open("garnier.txt", "r")

    helices = {}
    sheets = {}
    turns = {}
    coils = {}
    
    for line in f:
        if "Sequence:" in line:
            splitdata = line.split()
            seq = splitdata[2]
        elif "percent:" in line:
            percents = line.split()
            h = percents[3]
            e = percents[5]
            t = percents[7]
            c = percents[9]
            helices[seq] = h
            sheets[seq] = e
            turns[seq] = t
            coils[seq] = c

    return (helices, sheets, turns, coils)
            
def doMSA(finallistfile):
    args = ("t_coffee -in {0} -mode quickaln -output=score_ascii".format(finallistfile))
    os.system(args)
    
    f = open("sequences.score_ascii", "r")

    cons = {}
    
    for line in f:
        if "   :  " in line:
            splitdata = line.split()
            upid = splitdata[0]
            score = splitdata[2]
            cons[upid] = score
            
    return (cons)
      
    
def analyse(rxn, targ, pdir=0):
    
    rxnname = os.path.splitext(rxn)[0]
    csvname = rxnname.rsplit('/', 1)[-1]
    print ("Running quickRsim...")    
    
    (MnxSim, MnxDirPref, MnxDirUsed) = getMnxSim(rxn, pdir)
#    print(MnxSim)

    
    print ("Acquiring databases...")
    (sequence, names, descriptions, osource) = readFasta("data/seqs.fasta")
    
    with open('MnxToUprot.json') as f:
        MnxToUprot = json.load(f)
        
    with open ('upclst.json') as f2:
        upclst = json.load(f2)
        
    with open ('clstrep.json') as f3:
        clstrep= json.load(f3)

#    (clstup, upclst, clstrep) = clustar.readfile("cdhit_final.clstr")   # PUT THESE OUTSIDE, DO NOT REREAD EACH RUN
#    

#    MnxToUprot = readReacFile("reac_seqs.tsv")
    targplus = int(targ)*3
    list_mnx = sorted(MnxSim, key=MnxSim.__getitem__, reverse=True)[:int(targplus)]  #allow user to manipulate window of initial rxn id list
#    print (list_mnx)
    f = open("sequences.txt", 'w')
    print ("Creating initial MNX list...")
    targets = set()
#    global UprotToMnx
    UprotToMnx = {}
    
    # first creating fasta file, f, for further data extraction
    for x in list_mnx:
        up = MnxToUprot.get(x)  
        if up is not None:
            for y in up:
                if len(targets) >= int(targ):    # allow user to manipulate desired number of entries for resulting table
                    break
                else:
                    UprotToMnx[y] = x
                    try:
#                        cn = clustar.getClustNo(y)
#                        repid = clustar.getRepID(cn)
                        targets.add(y)
                    except KeyError:
                        pass 
  
    for t in targets:
        try: 
            seq = sequence[t]
            print ('>{0} \n{1}'.format(t, seq), file=f)
        except KeyError:
            pass

                    
    f.close()
    #analysis of FinalList of sequences
    (hydrop, weight, isoelec, polar) = pepstats("sequences.txt")
    (helices, sheets, turns, coils) = garnier("sequences.txt")
    cons = doMSA("sequences.txt")
    
    print ("Acquiring sequence properties...")
    # final table, do all data and value storing before this!
    rows = []
    for y in targets:
        try:
            desc = descriptions[y] 
            org = osource[y]
            mnx = UprotToMnx[y]
            cn = upclst.get(y)
            repid = clstrep[cn]
            rxnsimpre = float(MnxSim[mnx])
            rxnsim = float("{0:.5f}".format(rxnsimpre))
            conservation = float(cons[y])
#            placeholder = float(0)
            h = helices[y]
            e = sheets[y]
            t = turns[y]
            c = coils[y]
            w = weight[y]
            i = isoelec[y]
            pol = polar[y]
            
            rxndirused = MnxDirUsed[mnx]
            rxndirpref = MnxDirPref[mnx]
            
            rows.append( (y, desc, org, mnx, cn, repid, conservation, rxnsim, rxndirused, rxndirpref, h, e, t, c, w, i, pol) )
       
#            print ("{0}\t\t{1}\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{7}\t\t{8}\t\t{9}\t{10}\t\t{11}".format(repid, mnx, cn, rxndist, placeholder, h, e, t, c, w, i, pol))
        except KeyError:
            pass

    sortrows = sorted(rows, key = lambda x: (-x[7], -x[6]) )
    
#    f2 = open((os.path.join("uploads/table_"+rxnname+".txt")), 'w')
#    print ("Seq ID" + "\t\t" + "Rxn ID"+ "\t\t"+ "Clust. No." + "\t" + "Rep. ID."+ "\t" + "Rxn Sim."+ "\t" + "Consv. Score" + "\t" + "% helices" + "\t" + "% sheets" + "\t" + "% turns" + "\t\t" + "% coils" + "\t\t" + "Mol. Weight" + "\t" + "Isoelec. Point" + "\t" + "Polar %", file = f2)
#    for r in sortrows:
#        print ('\t'.join(map(str, r)))
#    f2.close()
    
    with open (os.path.join("results_"+csvname+".csv"), 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(('Seq. ID','Description', 'Organism Source', 'Rxn. ID', 'Clust. No.', 'Rep. ID.', 'Consv. Score', 'Rxn Sim.', "Direction Used", "Direction Preferred", '% helices', '% sheets', '% turns', '% coils', 'Mol. Weight', 'Isoelec. Point', 'Polar %'))
        for r in sortrows:
            writer.writerow(r)
    print ("CSV file created.")
    
    
def arguments():
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('rxn', 
                        help='Input rxn reaction file')
    parser.add_argument('-tar', type=float, default=20,
                        help='Number of targets to display in results [default = 20')
    parser.add_argument('-d', type=float, default=0,
                        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]')
    
    arg = parser.parse_args()
    return arg

if __name__ == '__main__':
    arg = arguments()
    
    analyse(arg.rxn, arg.tar, arg.d)













































