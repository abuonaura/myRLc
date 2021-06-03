import ROOT as r
import numpy as np
import os, sys

def CreateRootFiles():
    files2Hadd = ''
    for i in range(3):
        files2Hadd+=' prova_'+str(i)+'.root'
        f = r.TFile('prova_'+str(i)+'.root','RECREATE')
        t = r.TTree('t','t')
        entry = np.zeros(1, dtype=float)
        value = np.zeros(1, dtype=float)
        entry_b = t.Branch("entry",entry,"entry/D")
        value_b = t.Branch("value",value,"value/D")
        for j in range(100):
            entry[0] = i*100+j
            value[0] = r.TRandom3(0).Uniform(0,1)
            print(entry[0], value[0])
            t.Fill()
            #entry_b.Fill()
            #value_b.Fill()
        t.Write()
        f.Write()
        f.Close()
    return files2Hadd

def haddFiles(outfiles):
    outfile = 'provaHadded.root'
    os.system('hadd -f6 %s %s' %(outfile, outfiles))
    return outfile

def ReadHaddedFile(outfile):
    f = r.TFile(outfile,'READ')
    t = f.Get('t')
    print('Number of entries in the tree: ', t.GetEntries())
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        print(' Entry: ', i)
        print('  - Entry previous tree: ',t.entry)
        print('  - Value: ',t.value)



