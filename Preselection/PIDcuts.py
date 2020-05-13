import ROOT as r
from array import array

filedir = '$GANGAOUT/Datasets/'

def CheckPProbNNcut(tree,evt):
    tree.GetEntry(evt)
    if (tree.p_ProbNNp - tree.p_ProbNNk>0.):
        return 1
    else:
        return 0

def CheckMuPIDcut(tree,evt):
    tree.GetEntry(evt)
    if tree.mu_PIDmu>2 and tree.mu_PIDmu-tree.mu_PIDK>2 and tree.mu_PIDmu-tree.mu_PIDp>2.:
        return True
    else:
        return False

def GetPIDCutData(tree,evt,dtype):
    if dtype!='FakeMu' or dtype!='FakeMuSS':
        if CheckMuPIDcut(tree,evt):
            return 1
        else:
            return 0
    else:
        return 1

'''
def GetPIDCutMCfull(tree):
	if float(getattr(tree,'pi_PIDK_corr'))<2. and float(getattr(tree,'K_PIDK_corr'))>4. and float(getattr(tree,'p_PIDp_corr'))>0.:
		return True
	else:
		return False



def GetPIDCutMC(ifname):
	f = r.TFile(ifname[0:-5]+'_PIDCalib.root')
	t = f.Get('DecayTree')
'''

def CreatePIDcutTree(ifname,tree,dtype):
	of = r.TFile(ifname[0:-5]+'_PIDcut.root','RECREATE')
	ot = r.TTree('DecayTree','DecayTree')
	PIDcut = array('i',[0])
	ot.Branch('PIDcut',PIDcut,'PIDcut/I')
	if dtype!='MCfull' and dtype!='MCtrackeronly':
		for i in range(tree.GetEntries()):
			tree.GetEntry(i)
			if GetPIDCutData(tree,dtype):
				PIDcut[0]=1
			else:
				PIDcut[0]=0
			ot.Fill()
	of.Write()
	of.Close()
	return

if __name__ == "__main__":
	dtype = ['Data','DataSS','FakeMu','FakeMuSS','MCfull','MCtrackeronly']
	polarities=['MagUp','MagDown']
	for polarity in polarities:
		print(polarity)
		for dt in dtype:
			print('   ', dt)
			if dt!='MCfull' and dt!='MCtrackeronly':
				ifname = filedir+'Data/Lb_'+dt+'_'+polarity+'.root'
				infile = r.TFile(ifname,'READ')
				t = infile.Get('tupleout/DecayTree')
				CreatePIDcutTree(ifname,t,dt)
				infile.Close()


	


