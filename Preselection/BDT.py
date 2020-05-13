import ROOT as r
from array import array


myfiledir = '$GANGAOUT/Datasets/'
filedir_a = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/Data/'

def CreateBDTtreeData(dtype,polarity):
	f = r.TFile(filedir_a+dtype+'_'+polarity+'_2016_MVA.root','READ')
	t = f.Get('DecayTree')
	print(f, t)
	of = r.TFile(myfiledir+'Data/Lb_'+dtype+'_'+polarity+'_BDT.root', 'RECREATE')
	ot = r.TTree('DecayTree','DecayTree')
	bdt = array('f',[0.])
	LcMass = array('i',[0])
	ot.Branch('bdt',bdt,'bdt/F')
	ot.Branch('LcMass',LcMass,'LcMass/I')

	for i in range(t.GetEntries()):
		t.GetEntry(i)
		bdt[0], LcMass[0] = 0,0
		if t.Lc_M>2230 and t.Lc_M<2330:
			LcMass[0]=1
		bdt[0] = t.bdt
		#print(bdt[0])
		ot.Fill()
	of.Write()
	of.Close()
	f.Close()

if __name__ == "__main__":
	#dtype = ['Data','DataSS','FakeMu','FakeMuSS','MCfull','MCtrackeronly']
	dtype = ['Data','FakeMu','FakeMuSS']
	polarities=['MagUp','MagDown']
	for polarity in polarities:
		for dt in dtype:
			if dt!='MCfull' and dt!='MCtrackeronly':
				CreateBDTtreeData(dt,polarity)


	
