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



def readFasta(datadir, fileFasta):
    
    from Bio import SeqIO
    
    sequence = {}
    descriptions = {}
    osource = {}
    names = []
    seen = set()
    seen_add = seen.add
    
    for seq_record in SeqIO.parse(os.path.join(datadir, fileFasta), "fasta"):
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
    
def getMnxSim(rxn, p2env, datadir, outdir, drxn=0):
    cmd = [p2env, 'quickRsim.py', 
           os.path.join(datadir,'reac_prop.tsv'), os.path.join(datadir,'fp.npz'), '-rxn', rxn, '-out', os.path.join(outdir,'results_quickRsim.txt')]
    job = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = job.communicate()
#    print(out)
#    global mnx
#    mnx = []
    MnxSim = {}
    MnxDirPref = readRxnCons(os.path.join(datadir, "rxn_consensus_20160612.txt"))
    MnxDirUsed = {}
    
    if drxn==1:
        file = open(os.path.join(outdir, "results_quickRsim.txt"), 'r')
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

        file = open(os.path.join(outdir, "results_quickRsim.txt"), 'r')
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

   
def pepstats(file, outdir):
    outfile = os.path.join(outdir, "results.pepstats")
    args = ("pepstats -sequence {0} -outfile ".format(file) + outfile)
    os.system(args)
    
    f = open(outfile, "r")
    
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
        
def garnier(file, outdir):
    outfile = os.path.join(outdir, "garnier.txt")
    args = ("garnier -sequence {0} -outfile ".format(file) + outfile)
    os.system(args)
    
    f = open(outfile, "r")

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
            
def doMSA(finallistfile, outdir):
    outfile = os.path.join(outdir, "sequences.score_ascii")
    treefile = os.path.join(outdir, "sequences.dnd")
    args = ("t_coffee -in {0} -mode quickaln -output=score_ascii -outfile ".format(finallistfile) +outfile+ " -newtree "+treefile)
    os.system(args)
    
    f = open(outfile, "r")

    cons = {}
    
    for line in f:
        if "   :  " in line:
            splitdata = line.split()
            upid = splitdata[0]
            score = splitdata[2]
            cons[upid] = score
            
    return (cons)
      
    
def analyse(rxn, p2env, targ, datadir, outdir, csvfilename, pdir=0):
    
    datadir = os.path.join(datadir)
    outdir = os.path.join(outdir)
    
#    rxnname = os.path.splitext(rxn)[0]
#    csvname = rxnname.rsplit('/', 1)[-1]
    
    if csvfilename: 
        csvfilename = csvfilename
    else:
        csvfilename = "results_selenzy.csv"
        
    print ("Running quickRsim...")    
    (MnxSim, MnxDirPref, MnxDirUsed) = getMnxSim(rxn, p2env, datadir, outdir, pdir)
#    print(MnxSim)

    
    print ("Acquiring databases...")
    (sequence, names, descriptions, osource) = readFasta(datadir, "seqs.fasta")
    
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
    f = open(os.path.join(outdir, "sequences.txt"), 'w')
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
    (hydrop, weight, isoelec, polar) = pepstats(os.path.join(outdir, "sequences.txt"), outdir)
    (helices, sheets, turns, coils) = garnier(os.path.join(outdir, "sequences.txt"), outdir)
    cons = doMSA(os.path.join(outdir, "sequences.txt"), outdir)
    
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
    
    with open (os.path.join(outdir, csvfilename), 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(('Seq. ID','Description', 'Organism Source', 'Rxn. ID', 'Clust. No.', 'Rep. ID.', 'Consv. Score', 'Rxn Sim.', "Direction Used", "Direction Preferred", '% helices', '% sheets', '% turns', '% coils', 'Mol. Weight', 'Isoelec. Point', 'Polar %'))
        for r in sortrows:
            writer.writerow(r)
    print ("CSV file created.")
    
    
def arguments():
    parser = argparse.ArgumentParser(description='SeqFind script for Selenzy')
    parser.add_argument('rxn', 
                        help='Input rxn reaction file')
    parser.add_argument('p2env', 
                        help='Specify path to python 2 environment directory')
    parser.add_argument('-tar', type=float, default=20,
                        help='Number of targets to display in results [default = 20')
    parser.add_argument('-d', type=float, default=0,
                        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]')
    parser.add_argument('datadir',
                        help='specify data directory for required databases files, please end with slash')
    parser.add_argument('outdir',
                        help='specify output directory for all output files, including final CSV file, please end with slash')
    parser.add_argument('-outfile',
                        help='specify non-default name for CSV file output')
    
    arg = parser.parse_args()
    return arg

if __name__ == '__main__':
    arg = arguments()
    
    newpath = os.path.join(arg.outdir)
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    analyse(arg.rxn, arg.p2env, arg.tar, arg.datadir, arg.outdir, arg.outfile, arg.d)

#    from os import listdir
#    from os.path import isfile, join
#    
#    WORK = '/home/jerrywzy/Python Scripts/reactions'
#    
#    files = [f for f in listdir(WORK) if isfile(join(WORK, f))]
#
#    for f in files:
#        path = 'reactions/' + f
#        analyse(path, 100, 1)
#














