import ROOT as r
import preselection
import control_mode_JpsiK_fitnplot_Data
import JPSI

def apply_L0FilterAndPreselection( fileName_outOfGanga = "/disk/lhcb_data2/ibezshyi/DVntuple_DataTempor.root",treeName_outOfGanga = "TupleBToJpsi_K/DecayTree" ):
    file_no_selection = r.TFile(fileName_outOfGanga)
    tree_no_selection=file_no_selection.Get(treeName_outOfGanga)
    file_L0_preselection=r.TFile(fileName_outOfGanga[:-5]+"_L0_preselection.root","RECREATE")
    tree_L0_preselection = tree_no_selection.CopyTree(preselection.cut_string + " && (Bplus_L0MuonDecision_TOS || Bplus_L0Global_TIS) && Bplus_Hlt1TrackMuonDecision_TOS")
    tree_L0_preselection.Write()
    file_L0_preselection.Write()
    
def plotEfficiencies(efficiencies):
    #r.gStyle.SetOptStat(0)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionX().GetYaxis().SetRangeUser(2.1,2.8)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionY().GetYaxis().SetRangeUser(1.,2.9)
    efficiencies['Lb']['MC'][0]['TIS'].ProjectionX().GetYaxis().SetRangeUser(2.1,2.8)
    efficiencies['Lb']['MC'][0]['TIS'].ProjectionY().GetYaxis().SetRangeUser(1.,2.9)
    efficiencies['Lb']['MC'][1]['TIS'].SetLineColor(r.kBlue)
    efficiencies['JPSI']['MC'][0]['TISTOS'].SetLineColor(r.kGreen)
    efficiencies['JPSI']['MC'][1]['TISTOS'].SetLineColor(r.kGreen)
    efficiencies['JPSI']['MC'][0]['TIS'].SetLineColor(r.kOrange)
    efficiencies['JPSI']['MC'][1]['TIS'].SetLineColor(r.kOrange)
    efficiencies['Lb']['MC'][0]['TIS'].SetLineColor(r.kBlue)
    efficiencies['JPSI']['DATA'][0]['TISTOS'].SetLineColor(r.kViolet)
    efficiencies['JPSI']['DATA'][1]['TISTOS'].SetLineColor(r.kViolet)
    #efficiencies['JPSI']['MC'][0]['TISTOS'].SetMarkerColor(r.kGreen)
    #efficiencies['JPSI']['MC'][0]['TIS'].SetMarkerColor(r.kOrange)
    #efficiencies['Lb']['MC'][0]['TIS'].SetMarkerColor(r.kBlue)
    #efficiencies['JPSI']['DATA'][0]['TISTOS'].SetMarkerColor(r.kViolet)
    fN = r.TFile("efficiencies_JPSI.root","RECREATE")
          
    for iT in efficiencies.keys():
         for iD in efficiencies[iT].keys():
             for ik in xrange(len(efficiencies[iT][iD])):
                 for iM in efficiencies[iT][iD][ik].keys():
                     print iM, '  '
                     print efficiencies[iT][iD][ik][iM]
                     efficiencies[iT][iD][ik][iM].SetName(iT+"_"+iD+"_"+iM+"_"+str(ik))
                     efficiencies[iT][iD][ik][iM].Write()
    fN.Write()
    fN.Close()
    plots = {'cMethod':r.TCanvas('TISTOSvsTIS','TISTOS method vs. TIS'),'cCompat':r.TCanvas('JPSIvsLb','TIS for JPSI vs. Lb'),'cSim':r.TCanvas('DATAvsMC','TISTOS DATA vs. MC')}
    for iC in plots.keys(): plots[iC].Divide(2,2)
    #legend =r.TLegend(0.1,0.7,0.48,0.9)
    #legend.AddEntry(efficiencies['JPSI']['DATA'][0]['TISTOS'],"Data: TISTOS/TOS (J/\Psi)","lep")
    #legend.AddEntry(efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionX(),"MC: TISTOS/TOS (J/\Psi)","lep")
    #legend.AddEntry(efficiencies['JPSI']['MC'][0]['TIS'],"TIS/ALL (J/\Psi)","lep")
    #legend.AddEntry(efficiencies['Lb']['MC'][0]['TIS'],"TIS/ALL (L_{c}\mu)","lep")
    
    plots['cMethod'].cd(1)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionX().Draw()
    efficiencies['JPSI']['MC'][0]['TIS'].ProjectionX().Draw('same')
    plots['cMethod'].cd(2)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionY().Draw()
    efficiencies['JPSI']['MC'][0]['TIS'].ProjectionY().Draw('same')
    plots['cMethod'].cd(3)
    efficiencies['JPSI']['MC'][1]['TISTOS'].Draw()
    efficiencies['JPSI']['MC'][1]['TIS'].Draw('same')
  #  plots['cMethod'].cd(4)
 #   legend.Draw()
    
    
    plots['cCompat'].cd(1)
    efficiencies['JPSI']['MC'][0]['TIS'].ProjectionX().Draw()
    efficiencies['Lb']['MC'][0]['TIS'].ProjectionX().Draw('same')
    plots['cCompat'].cd(2)
    efficiencies['JPSI']['MC'][0]['TIS'].ProjectionY().Draw()
    efficiencies['Lb']['MC'][0]['TIS'].ProjectionY().Draw('same')
    plots['cCompat'].cd(3)
    efficiencies['JPSI']['MC'][1]['TIS'].Draw()
    efficiencies['Lb']['MC'][1]['TIS'].Draw('same')
    #plots['cCompat'].cd(4)
    #legend.Draw()
        
    plots['cSim'].cd(1)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionX().Draw()
    efficiencies['JPSI']['DATA'][0]['TISTOS'].ProjectionX().Draw('same')
    plots['cSim'].cd(2)
    efficiencies['JPSI']['MC'][0]['TISTOS'].ProjectionY().Draw()
    efficiencies['JPSI']['DATA'][0]['TISTOS'].ProjectionY().Draw('same')
    plots['cSim'].cd(3)
    efficiencies['JPSI']['MC'][1]['TISTOS'].Draw()
    efficiencies['JPSI']['DATA'][1]['TISTOS'].Draw('same')
#    plots['cSim'].cd(4)
  #  legend.Draw()
    for iC in plots.keys(): 
        plots[iC].SaveAs(iC+".png")
        plots[iC].SaveAs(iC+".root")


        
    
    
def main(checks=1):
    MC_fileName = "/disk/lhcb_data2/ibezshyi/DVntuple_MC.root"
    Data_fileName = "/disk/lhcb_data2/ibezshyi/DVntuple_DataTempor.root"
    Lb_fileName = "/disk/lhcb_data2/ibezshyi/Lb_Lcmunu_MagUp.root"
    
# we apply preselection based on L0 trigger and basic PID and mass decisions for MC...
    apply_L0FilterAndPreselection(fileName_outOfGanga = MC_fileName)
# ... and data
    apply_L0FilterAndPreselection(fileName_outOfGanga = Data_fileName,treeName_outOfGanga = "DecayTree")
# we calculate s_weights for our tuple
    control_mode_JpsiK_fitnplot_Data.main(Data_fileName[:-5]+"_L0_preselection.root")
# we repeat for MC just with a reason to compare a mass fits
    control_mode_JpsiK_fitnplot_Data.main(MC_fileName[:-5]+"_L0_preselection.root")
# we create P_t/P_z efficiency histogram
    efficiencies = {'JPSI':{'MC':[],'DATA':[]},'Lb':{'MC':[],'DATA':[]}}
    efficiencies['JPSI']['DATA'] = JPSI.TISTOSvsREAL(Data_fileName[:-5]+"_L0_preselection.root",'DecayTree', dataormc = 'sweight')
    if checks:
        # JPSI mc
        efficiencies['JPSI']['MC'] = JPSI.TISTOSvsREAL(MC_fileName[:-5]+"_L0_preselection.root",'DecayTree')
        # Lb mc
        efficiencies['Lb']['MC'] = JPSI.TISTOSvsREAL(Lb_fileName,'tupleout/DecayTree','Lb','Lc_L0HadronDecision_TOS')
        plotEfficiencies(efficiencies)
        
