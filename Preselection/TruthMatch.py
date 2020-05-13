import ROOT as r
from array import array

def TruthMatch(tree):
    if int(getattr(tree,'Lc_BKGCAT'))<30 and int(getattr(tree,'Lb_BKGCAT'))<50:
        flag = True
    else:
        flag = False
    return flag

def TruthMatchLambdab(tree):
    if abs(int(getattr(tree,'Lb_TRUEID')))==5122:
        flag=True
    else:
        flag=False
    return flag

def TruthMatchLambdac(tree):
    if abs(int(getattr(tree,'Lc_TRUEID')))==4122:
        flag=True
    else:
        flag=False
    return flag

def TruthMatchMuonPdg(tree):
    if abs(int(getattr(tree,'mu_TRUEID')))==13:
        flag=True
    else:
        flag = False
    return flag

def TruthMatchMuonCharm(tree):
    if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==431 or abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==421 or  abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==411:
        flag = True
    else:
        flag=False
    return flag

def TruthMatchMuonLb(tree):
    if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==5122:
        flag = True
    else:
        flag=False
    return flag

def TruthMatchMuonTau(tree):
    if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==15:
        flag = True
    else:
        flag=False
    return flag

def GetTruthMatchCut(sample):
    if sample in ['Lb_LcDs','Lb_Lc2593Ds','Lb_Lc2625Ds']:
        cut = r.TCut('Lc_BKGCAT<30&&Lb_BKGCAT<50&&abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122&&abs(mu_TRUEID)==13&&(abs(mu_MC_MOTHER_ID)==431||abs(mu_MC_MOTHER_ID)==421||abs(mu_MC_MOTHER_ID)==411)')
    elif sample in ['Lb_Lctaunu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
        cut = r.TCut('Lc_BKGCAT<30&&Lb_BKGCAT<50&&abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122&&abs(mu_TRUEID)==13&&abs(mu_MC_MOTHER_ID)==15')
    else:
        cut = r.TCut('Lc_BKGCAT<30&&Lb_BKGCAT<50&&abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122&&abs(mu_TRUEID)==13&&abs(mu_MC_MOTHER_ID)==5122')
    return cut

def GetFinalTruthMatch(sample,tree):
	if TruthMatch(tree) and TruthMatchLambdab(tree) and TruthMatchLambdac(tree) and TruthMatchMuonPdg(tree):
		if sample in ['Lb_LcDs','Lb_Lc2593Ds','Lb_Lc2625Ds']:
			if TruthMatchMuonCharm(tree):
				return True
			else:
				return False
		elif sample in ['Lb_Lctaunu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
			if TruthMatchMuonTau(tree):
				return True
			else:
				return False
		else:
			if TruthMatchMuonLb(tree):
				return True
			else:
				return False


def CreateTruthMatchTree(ifname,tree,sample):
	print(tree)
	of = r.TFile(ifname[0:-5]+'_TrutMatch.root','RECREATE')
	ot = r.TTree('DecayTree','DecayTree')
	TM_comb = array( 'i', [ 0 ] )
	TM_Lb = array( 'i', [ 0 ] )
	TM_Lc = array( 'i', [ 0 ] )
	TM_Mu = array( 'i', [ 0 ] )
	TM_MuCharm = array( 'i', [ 0 ] )
	TM_MuTau = array( 'i', [ 0 ] )
	TM_MuLb = array( 'i', [ 0 ] )
	TM_Final = array( 'i', [ 0 ] )
	ot.Branch( 'TruthMatch_comb', TM_comb, 'TruthMatch_comb/I' )
	ot.Branch( 'TruthMatch_Lambdab', TM_Lb, 'TruthMatch_Lambdab/I' )
	ot.Branch( 'TruthMatch_Lambdac', TM_Lc, 'TruthMatch_Lambdac/I' )
	ot.Branch( 'TruthMatch_Mu', TM_Mu, 'TruthMatch_Mu/I' )
	ot.Branch( 'TruthMatch_MuFromCharm', TM_MuCharm,'TruthMatch_MuFromCharm/I')
	ot.Branch( 'TruthMatch_MuFromTau', TM_MuTau,'TruthMatch_MuFromTau/I')
	ot.Branch( 'TruthMatch_MuFromLb', TM_MuLb,'TruthMatch_MuFromLb/I')
	ot.Branch( 'TruthMatch_Final', TM_Final,'TruthMatch_Final/I')

	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		TM_comb[0]=0
		TM_Lb[0]=0
		TM_Lc[0]=0
		TM_Mu[0]=0
		TM_MuCharm[0]=0
		TM_MuTau[0]=0
		TM_MuLb[0]=0
		TM_Final[0]=0
		if TruthMatch(tree) and TruthMatchLambdab(tree) and TruthMatchLambdac(tree) and TruthMatchMuonPdg(tree):
			TM_comb[0]=1
			TM_Lb[0]=1
			TM_Lc[0]=1
			TM_Mu[0]=1
			if sample in ['Lb_LcDs','Lb_Lc2593Ds','Lb_Lc2625Ds']:
				if TruthMatchMuonCharm(tree):
					TM_MuCharm[0]=1
			elif sample in ['Lb_Lctaunu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
				if TruthMatchMuonTau(tree):
					TM_MuTau[0]=1
			else:
				if TruthMatchMuonLb(tree):
					TM_MuLb[0]=1
		if GetFinalTruthMatch(sample,tree):
			TM_Final[0]=1
		ot.Fill()

	of.Write()
	of.Close()

if __name__ == "__main__":
	filedir = '/disk/users/iaroslava/'
	ifname = filedir+'Lc_MVAupdate_MagUp.root'
	f = r.TFile(ifname,'READ')
	t = f.Get('tupleout/DecayTree')
	print(f, t)
	CreateTruthMatchTree(ifname, t,'Lb_Lcmunu')
