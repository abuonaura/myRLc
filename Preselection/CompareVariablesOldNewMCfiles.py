import ROOT as r
import os

directory = {'old':'/disk/lhcb_data2/RLcMuonic2016/MC_full_trueTrigger/','new':'/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/'}
mcsamples  = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds']
polarities = ['MagUp','MagDown']

def PutTogetherPolarityHistos(h):
    h_new = {}
    for mcsample in mcsamples:
        h_new[mcsample] = h[mcsample]['MagUp']
        h_new[mcsample].Add(h[mcsample]['MagDown'])
        h_new[mcsample].SetDirectory(0)
    return h_new

def GetHistograms(variable):
    hold = {mcsample: {polarity: r.TH1F('hold_'+variable+'_'+mcsample+'_'+polarity,'',40,-2,12) for polarity in polarities} for mcsample in mcsamples}
    hnew = {mcsample: {polarity: r.TH1F('hnew_'+variable+'_'+mcsample+'_'+polarity,'',40,-2,12) for polarity in polarities} for mcsample in mcsamples}
    for mcsample in mcsamples:
        for polarity in polarities:
            fold = r.TFile(directory['old']+mcsample+'_'+polarity+'_full.root','READ')
            told = fold.Get('tupleout/DecayTree')
            name = mcsample+'_'+polarity
            told.Draw(variable)
            hold[mcsample][polarity] = r.gPad.GetPrimitive('hist')
            fnew = r.TFile(directory['new']+mcsample+'_'+polarity+'_full.root','READ')
            tnew = fnew.Get('tupleout/DecayTree')
            tnew.Draw(variable)
            hnew[mcsample][polarity] = r.gPad.GetPrimitive('hist')
            hold[mcsample][polarity].SetDirectory(0)
            hnew[mcsample][polarity].SetDirectory(0)
    hold1 = PutTogetherPolarityHistos(hold)
    hnew1 = PutTogetherPolarityHistos(hnew)
    return hold1,hnew1

def PlotSuperimposed(hold, hnew):
    c = r.TCanvas()
    hold.SetLineColor(r.Black)
    hnew.SetLineColor(r.kOrange+2)
    hold.Draw('hist')
    hnew.Draw('hist sames')
    l = r.TLegend()
    l.AddEntry(hold,'old file','l')
    l.AddEntry(hnew,'new file','l')
    l.Draw('same')
    return c

def CompareAllVariables():
    parentdir = os.getcwd()
    mcsamples = ['Lb_LcDs']
    for mcsample in mcsamples:
        for polarity in polarities:
            plotdir = 'ComparisonPlots/'+mcsample+'/'+polarity
            path = os.path.join(parentdir, plotdir)
            try:
                os.makedirs(path)
            except OSError:
                print('Directory already existing')
            else:
                print("Successfully created the directory %s" % path)
            fold = r.TFile(directory['old']+mcsample+'_'+polarity+'_full.root','READ')
            told = fold.Get('tupleout/DecayTree')
            fnew = r.TFile(directory['new']+mcsample+'_'+polarity+'_full.root','READ')
            tnew = fnew.Get('tupleout/DecayTree')

            branches = told.GetListOfBranches()
            for i in range(branches.GetEntries()):
                cname = 'c_'+branches.At(i).GetName()+'_'+mcsample+'_'+polarity 
                c = r.TCanvas(cname,cname,500,500)
                told.Draw(branches.At(i).GetName())
                tnew.Draw(branches.At(i).GetName(),'','sames')
                r.gPad.GetListOfPrimitives().At(1).SetLineColor(2);
                c.SaveAs(path+'/'+cname+'.png')

CompareAllVariables()
            
'''
told.Draw(variable)
tnew.Draw(variable)
'''


        

