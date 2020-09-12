import ROOT as r


filedir = '/disk/lhcb_data2/RLcMuonic2016/DataPreselectedFull/'
datatype = ['Data','DataSS','FakeMu','FakeMuSS']
polarities = ['MagUp','MagDown']

for dt in datatype:
    nsel, nsel_sw, nevt, nevt_sw=0,0,0,0
    nsel_ai, nsel_ai_sw=0,0
    for polarity in polarities:
        fname = 'Lb_'+dt+'_'+polarity+'_preselected_full_sw.root'
        f = r.TFile(filedir+fname,'Read')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            nevt+=1
            nevt_sw+=t.sw_sig
            if t.Lb_ISOLATION_BDT<0.35:
                nsel+=1
                nsel_sw+=t.sw_sig
            if t.Lb_ISOLATION_BDT>0.35 and t.Lb_ISOLATION_BDT2>0.2:
                nsel_ai+=1
                nsel_ai_sw+=t.sw_sig

    print('---- Datatype: ',dt)
    print('Total number of events: ',nevt)
    print('Requiring IsolationBDT<0.35: ' )
    print('Number of selected evts: '+str(nsel)+'  eff: {:.2f} %'.format(nsel*100./nevt))
    print('(With sweights) Number of selected evts: '+str(nsel_sw)+'  eff: {:.2f} %'.format(nsel_sw*100./nevt_sw))
    print('Requiring IsolationBDT>0.35 && IsolationBDT2>0.2: ' )
    print('Number of selected evts: '+str(nsel_ai)+'  eff: {:.2f} %'.format(nsel_ai*100./nevt))
    print('(With sweights) Number of selected evts: '+str(nsel_ai_sw)+'  eff: {:.2f} %'.format(nsel_ai_sw*100./nevt_sw))


