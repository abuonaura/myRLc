import ROOT as r


datatypes = ['Data','FakeMu','DataSS','FakeMuSS']
#datatypes = ['Data']
polarities=['MagUp','MagDown']
suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}
suffix_sw = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}


def retrieveSignalYield(fname):
    f = r.TFile(fname,'READ')
    t = f.Get('DecayTree')
    h = r.TH1F('h','',100,2230,2330)
    t.Draw('Lc_M>>h','sw_sig*(Lb_M>0)')
    h = r.gPad.GetPrimitive('h')
    yd = h.Integral()
    return yd

def retrieveBkgYield(fname):
    f = r.TFile(fname,'READ')
    t = f.Get('DecayTree')
    h = r.TH1F('h','',100,2230,2330)
    t.Draw('Lc_M>>h','sw_bkg*(Lb_M>0)')
    h = r.gPad.GetPrimitive('h')
    yd = h.Integral()
    return yd

def PrintResults(yields):
    print('        Signal:')
    for sample in suffix:
        print('        - Sample: ', sample) 
        print('                MagUp: {:.2f},  MagDown: {:.2f}'.format(yields['Data']['MagUp'][sample], yields['Data']['MagDown'][sample]))
    print('        FakeMu:')
    for sample in suffix:
        print('        - Sample: ', sample) 
        print('                MagUp: {:.2f},  MagDown: {:.2f}'.format(yields['FakeMu']['MagUp'][sample], yields['FakeMu']['MagDown'][sample]))
    print('        Same Sign:')
    for sample in suffix:
        print('        - Sample: ', sample) 
        print('                MagUp: {:.2f},  MagDown: {:.2f}'.format(yields['DataSS']['MagUp'][sample], yields['DataSS']['MagDown'][sample]))
    print('        FakeMu Same Sign:')
    for sample in suffix:
        print('        - Sample: ', sample) 
        print('                MagUp: {:.2f},  MagDown: {:.2f}'.format(yields['FakeMuSS']['MagUp'][sample], yields['FakeMuSS']['MagDown'][sample]))



if __name__=='__main__':
    s_yields = {}
    b_yields = {}
    for datat in datatypes:
        s_yields[datat] ={}
        b_yields[datat] ={}
        for polarity in polarities:
            s_yields[datat][polarity] = {}
            b_yields[datat][polarity] = {}
            for sample in suffix_sw:
                #print(suffix_sw[sample])
                if sample == 'full' or sample == 'iso':
                    datadir = '$FILEDIR/Data/'
                else:
                    datadir = '$FILEDIR/ControlSamples/'
                fname =datadir+'Lb_'+datat+'_'+polarity+'_reduced_preselected'+suffix_sw[sample]
                s_yields[datat][polarity][sample]=retrieveSignalYield(fname)
                b_yields[datat][polarity][sample]=retrieveBkgYield(fname)
    print('----------   SIGNAL YIELDS   ------------')
    PrintResults(s_yields)
    print('----------   BKG YIELDS   ------------')
    PrintResults(b_yields)



