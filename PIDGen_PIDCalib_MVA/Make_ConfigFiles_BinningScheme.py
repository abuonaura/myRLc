import yaml
import os

UraniaDir = sys.argv[1]
#UraniaDir = '/home/hep/amathad/Packages/UraniaDev_v8r0/' 

cut          = {}
cut["Mu"]    = "Brunel_P>3000&&nSPDHits<600&&Brunel_IPCHI2>16&&Brunel_TRACK_GHOSTPROB<0.5"
cut["Pi"]    = "Brunel_PT>300&&Brunel_P>2000&&nSPDHits<600&&Brunel_IPCHI2>9&&Brunel_TRACK_GHOSTPROB<0.5"
cut["K"]     = "Brunel_PT>300&&Brunel_P>2000&&nSPDHits<600&&Brunel_IPCHI2>9&&Brunel_TRACK_GHOSTPROB<0.5"
cut["P"]     = "Brunel_PT>300&&Brunel_P>2000&&nSPDHits<600&&Brunel_IPCHI2>9&&Brunel_TRACK_GHOSTPROB<0.5"

pidcut       = {}
pidcut["Mu"] = "DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1"
pidcut["Pi"] = "DLLK<2"
pidcut["K"]  = "DLLK>4"
pidcut["P"]  = "DLLp>0"

varminmax = {}
varminmax["Brunel_P"]       = (0, 200000, 100)
varminmax["Brunel_PT"]      = (0,  60000, 100)
varminmax["nTracks_Brunel"]        = (0,  700  ,  5)

for data in ["Turbo16"]:
    for magtype in ["MagUp"]:
        for particle in ["Pi", "P", "K", "Mu"]:
            for varname in ["Brunel_P", "Brunel_PT", "nTracks_Brunel"]:
                config = {}
                config["sampleVersion"             ] = data
                config["magnetPolarity"            ] = magtype
                config["particleName"              ] = particle
                config["priorCut"                  ] = cut[particle]
                config["pidCut"                    ] = pidcut[particle]
                config["varName"                   ] = varname
                config["outputFile"                ] = "binoutput/binning-"+data+".py"
                config["minimum"                   ] = varminmax[varname][0]
                config["maximum"                   ] = varminmax[varname][1]
                config["minimumBinWidth"           ] = varminmax[varname][2]
                config["delta"                     ] = 1
                config["nSigma"                    ] = 5
                config["schemeName"                ] = 'binning-'+particle+'-'+data+'-'+magtype
                config["numberOfInitialBins"       ] = 100 
                config["startWithIsopopulatedBins" ] = True
                config["minRun"                    ] = None
                config["maxRun"                    ] = None
                config["maxFiles"                  ] = 100 #used
                config["mergeBelow"                ] = -1000000000
                config["mergeAbove"                ] = 1000000000
                #print(config)
                if not os.path.isdir('config_files'): os.mkdir('config_files')
                with open("config_files/config-"+particle+"-"+data+"-"+magtype+"-"+varname+".yml", "w") as outfile: yaml.dump(config, outfile, default_flow_style=False)

                command = "bash "+UraniaDir+"run python "+UraniaDir+"PIDCalib/PIDPerfScripts/scripts/python/BinningOptimizer/binningPID.py "+"config_files/config-"+particle+"-"+data+"-"+magtype+"-"+varname+".yml"
                print(command)
                os.system(command)
