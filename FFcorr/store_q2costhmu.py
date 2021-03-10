#!/bin/python

import numpy as np
#import root_pandas as rpd
import os

def LorentzBoost(vector, boostvector):
    boost = SpatialComponents(boostvector)
    b2    = ScalarProduct(boost, boost)
    gamma = 1./np.sqrt(1.-b2)
    gamma2= (gamma-1.0)/b2
    ve    = TimeComponent(vector)
    vp    = SpatialComponents(vector)
    bp    = ScalarProduct(vp, boost)
    vp2   = vp + Scalar(gamma2*bp + gamma*ve)*boost
    ve2   = gamma*(ve + bp)
    return LorentzVector(vp2, ve2)

def BoostToRest(vector, boostvector):
    boost = -SpatialComponents(boostvector)/Scalar(TimeComponent(boostvector))
    return LorentzBoost(vector, boost)

def UnitVector(vec):
    return vec/Scalar(Norm(vec))

def Norm(vec):
    return np.sqrt(np.sum(vec*vec, 1))

def MetricTensor():
    return np.array([-1., -1., -1., 1.])

def Mass(vector):
    return np.sqrt(np.sum(vector*vector*MetricTensor(), 1))

def ScalarProduct(vec1, vec2):
    return np.sum(vec1*vec2, 1)

def Scalar(x): 
    return np.stack([x], axis=1)

def Vector(x, y, z): 
    return np.stack([x, y, z], axis=1)

def LorentzVector(space, time): 
    return np.concatenate([space, np.stack([time], axis=1)], axis=1)

def TimeComponent(vector): 
    return vector[:, 3]

def SpatialComponents(vector): 
    return vector[:, 0:3]

def return_phasespace(PLb_lab, PLc_lab, PLepton_lab):
    PW_lab       = PLb_lab - PLc_lab
    PLb_Wlab     = BoostToRest(PLb_lab    , PW_lab) #Boost Lb to W rest frame
    PLepton_Wlab = BoostToRest(PLepton_lab, PW_lab) #Boost Lepton to W rest frame
    q2           = Mass(PW_lab)**2  #evaluate q^2
    #Evaluate costhl as the angle b/w opposite direction of Lb in W rest frame and Lepton in W rest frame
    costhl       = ScalarProduct(UnitVector(SpatialComponents(PLepton_Wlab)), -1. * UnitVector(SpatialComponents(PLb_Wlab))) 
    return np.nan_to_num(q2, nan=0.), np.nan_to_num(costhl, nan=10.)

def ComputeQ2andCosTheta():
    fname      = ['/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_Lcmunu_MagUp_full.root']
    treename   = 'tupleout/DecayTree'
    columns    = ['*TRUE*']
    #do not need cut
    cut        = 'abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122&&abs(p_TRUEID)==2212&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(mu_TRUEID)==13&&Lb_BKGCAT<60'
    parentname = 'Lb'
    dauglc     = 'Lc' #for Lc* samples this name should be replaced by Lc2595, Lc2625 etc.
    dauglep    = 'mu' #for tau samples this name should be replaced by tau
    if os.path.exists('./Lb_Lcmunu_MagUp_full_new.root'):
        print('Removing file ./Lb_Lcmunu_MagUp_full_new.root')
        os.remove('./Lb_Lcmunu_MagUp_full_new.root')
    
    count = 0
    for df in rpd.read_root(fname, columns = columns,key=treename,where=cut, chunksize=50000):
        print(df.shape)
        pxlb      = df[parentname+'_TRUEP_X']; pxlc = df[dauglc+'_TRUEP_X']; pxl = df[dauglep+'_TRUEP_X']; pxnu = pxlb - pxlc - pxl
        pylb      = df[parentname+'_TRUEP_Y']; pylc = df[dauglc+'_TRUEP_Y']; pyl = df[dauglep+'_TRUEP_Y']; pynu = pylb - pylc - pyl
        pzlb      = df[parentname+'_TRUEP_Z']; pzlc = df[dauglc+'_TRUEP_Z']; pzl = df[dauglep+'_TRUEP_Z']; pznu = pzlb - pzlc - pzl
        pelb      = df[parentname+'_TRUEP_E']; pelc = df[dauglc+'_TRUEP_E']; pel = df[dauglep+'_TRUEP_E']; penu = pelb - pelc - pel
        PLc_lab   = LorentzVector(Vector(pxlc, pylc, pzlc), pelc) #Format of LorentzVector(Vector(X,Y,Z), E)
        Pl_lab    = LorentzVector(Vector(pxl , pyl , pzl ), pel )
        PNu_lab   = LorentzVector(Vector(pxnu, pynu, pznu), penu)
        PLb_lab   = PLc_lab + Pl_lab + PNu_lab
        qsq, cthl = return_phasespace(PLb_lab, PLc_lab, Pl_lab)
        df['Q2TRUE']     = qsq
        df['COSTHLTRUE'] = cthl
        df.to_root('Lb_Lcmunu_MagUp_full_new.root', key='DecayTree', mode='a')
        #count += 1
        #if count == 10: exit(1)
        return
