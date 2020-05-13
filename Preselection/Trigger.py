import ROOT as r
from array import array


def L0TriggerData(t,evt):
	t.GetEntry(evt)
	if t.Lb_L0Global_TIS==1 or t.Lb_L0HadronDecision_TOS==1:
		return 1
	else:
		return 0

def HLT1TriggerData(t,evt):
	t.GetEntry(evt)
	if t.Lc_Hlt1TrackMVADecision_TOS==1 or t.Lc_Hlt1TwoTrackMVADecision_TOS==1:
		return 1
	else:
		return 0

def HLT2TriggerData(t,evt,dtype):
	t.GetEntry(evt)
	if dtype =='Data' or dtype == 'DataSS':
		if t.Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==1:
			return 1
		else:
			return 0
	if dtype =='FakeMu' or dtype=='FakeMuSS':
		if t.Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS==1:
			return 1
		else:
			return 0

def GetTriggerData(L0, HLT1, HLT2):
	if L0==1 and HLT1==1 and HLT2==1:
		return 1
	else:
		return 0



def CreateTriggerTreeData(ifname,tree,dtype):
	of = r.TFile(ifname[0:-5]+'_Trigger.root','RECREATE')
	ot = r.TTree('DecayTree','DecayTree')
	L0 = array( 'i', [ 0 ] )
	HLT1 = array( 'i', [ 0 ] )
	HLT2 = array( 'i', [ 0 ] )
	TFinal = array('i',[0])
	ot.Branch('L0',L0,'L0/I')
	ot.Branch('HLT1',HLT1,'HLT1/I')
	ot.Branch('HLT2',HLT2,'HLT2/I')
	ot.Branch('Trigger_Final',TFinal,'Trigger_Final/I')

	tree.SetBranchStatus('*',0)
	tree.SetBranchStatus('Lb_L0Global_TIS',1)
	tree.SetBranchStatus('Lb_L0HadronDecision_TOS',1)
	tree.SetBranchStatus('Lc_Hlt1TrackMVADecision_TOS',1)
	tree.SetBranchStatus('Lc_Hlt1TwoTrackMVADecision_TOS',1)
	if dtype =='Data' or dtype == 'DataSS':
		tree.SetBranchStatus('Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS',1)
	else:
		tree.SetBranchStatus('Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS',1)
	
	for i in range(tree.GetEntries()):
		L0[0], HLT1[0], HLT2[0], TFinal[0] = 0,0,0,0
		if L0Trigger(tree,i): L0[0]=1
		if HLT1Trigger(tree,i): HLT1[0]=1
		if HLT2Trigger(tree,i,dtype): HLT2[0]=1
		if TriggerFinal(L0[0],HLT1[0],HTL2[0]):TFinal[0]=1
		ot.Fill()
	of.Write()
	of.Close()


if __name__ == "__main__":
	filedir = '$GANGAOUT/Datasets/'
	dtype = ['Data','DataSS','FakeMu','FakeMuSS']
	polarities=['MagUp','MagDown']
	for polarity in polarities:
		for dt in dtype:
			if dt!='MCfull' and dt!='MCtrackeronly':
				ifname = filedir+'Data/Lb_'+dt+'_'+polarity+'.root'
				f = r.TFile(ifname,'READ')
				t = f.Get('tupleout/DecayTree')
				print(f, t)
				CreateTriggerTreeData(ifname,t,dt)



