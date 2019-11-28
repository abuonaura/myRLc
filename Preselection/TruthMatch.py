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


