# RelInfoHLT1Emulation

This tool can be used to save all the informations needed to emulate the HLT1 response offline.
For the moment it contains all the variables needed in order to emulate the following HLT1 lines:
    - HLT1TrackMVA
    - HLT1TwoTrackMVA

## Compilation

### Create a new Developement DaVinci Project
For the moment I have tested this tool with the version of DaVinci we are using in the RDplus analysis, that is `v42r6p1`

```
lb-dev DaVinci/v42r6p1
cd DaVinciDev_v42r6p1
```

#### Check out the packages you would like to modify
```
git lb-use Phys
git lb-checkout Phys/v23r6p1 Phys/LoKiPhys
git lb-checkout Phys/v23r6p1 Phys/DaVinciTypes 
git lb-checkout Phys/v23r6p1 Phys/RelatedInfoTools
```

### Modify the packages
```
cp ../RelInfoHLT1Emulation.* Phys/RelatedInfoTools/src/
cp ../RelatedInfoNamed.h     Phys/DaVinciTypes/Kernel/
```

### compile
```
make 
```

## Usage
The variables that are evaluated with this tuple will be available in DaVinci with the use of the `RELINFO` LoKi functor. The variables are the following:

  - HEAD_NDAUGHTERS             : Number of stable aughters found in the decay. maximum saved is MAXNUMBER
  - HEAD_NCOMBINATIONS          : Number of ombinations. Maximum is [MAXNUMBER*(MAXNUMBER-1)]/2
  - PT_DAU_i               : PT of the i^th stable aughter
  - P_DAU_i                : P of the i^th stable aughter
  - IPCHI2_OWNPV_DAU_i     : IPChi2 of the i^th aughter wrt its best PV
  - TRACK_GHOSTPROB_DAU_i  : Ghost Probability of he i^th daughter
  - TRACK_CHI2_DAU_i       : Track Chi2 of the ^th daughter
  - TRACK_NDOF_DAU_i       : Track NDoF of the ^th daughter
  - VERTEX_CHI2_COMB_i_j   : Vertex Chi2 of the -j combination
  - VERTEX_NDOF_COMB_i_j   : Vertex NDoF of the -j combination
  - ETA_COMB_i_j           : Pseudorapidity of the -j combination
  - MCORR_OWNPV_COMB_i_j   : Corrected Mass of the -j combination using the OWNPV
  - SUMPT_COMB_i_j         : Sum of the daughters' T of the i-j combination 
  - DIRA_OWNPV_COMB_i_j    : DIRA of the i-j ombination using the OWNPV
  - DOCA_COMB_i_j          : DOCAChi2  of the aughters of the i-j combination
  - VDCHI2_OWNPV_COMB_i_j  : VDCHI2 of the i-j ombination wrt the OWNPV 
  - IPCHI2_ONWPV_COMB_i_j  : IPCHI2 of the i-j ombination wrt the OWNPV
  - NLT_OWNPV_COMB_i_j     : Nlt variable used in the HLT1, wrt the OWNPV

The  NLT_OWNPV of a particle is the number of its daughters which have an IP_\chi^2 less than a given value. If not configured, this value will be 16.

In order to use this tool you need to firstly import it into DaVinci and set it up.

```
from Configurables import RelInfoHLT1Emulation

relinfo_HLT1.addTool(RelInfoHLT1Emulation, "RelInfoHLT1Emulation")
relinfo_HLT1.Tool = "RelInfoHLT1Emulation"
relinfo_HLT1.Location = 'HLT1Emulation'
relinfo_HLT1.Inputs = [selection_output]
hlt1_emu = getattr(relinfo_HLT1, "RelInfoHLT1Emulation")
hlt1_emu.nltValue = 16.

DaVinci().appendToMainSequence([relinfo_HLT1])
```

The `Inputs` member should point to the TES locations of the particles to which you want
to apply the RelatedInfo Tool. 

The `nltValue` members is the cut-value with which the NLT variable is evaluated. If not set this will be set to be 16.

The values can then be read back with a `LoKi::Hybrid::TupleTool`. The values are saved in a TES location that is the same for the Particles TES you are evaluating the Tool on, with `Particles` substituted with `HLT1Emulation`. For example, assume that tuple is your Tuple, and you want to save some of the variables into it:

```
relinfo_output = selection_output.replace("Particles", "HLT1Emulation")

tt_HLT1Emulation = tuple.addTupleTool("LoKi::Hybrid::TupleTool/Hlt1Emulation")
tt_HLT1Emulation.Preambulo = []
tt_HLT1Emultion.Variables = {
  "NDAUGHTERS" : "RELINFO('"+relinfo_output+"', 'NDAUGHTERS', 0)",
  "PT_DAU_1"   : "RELINFO('"+relinfo_output+"', 'PT_DAU_1'  , 0)"
}
```

Next you might want to evaluate some MVA on the output of this tool directly in DaVinci. This can be done with the `MVADictHelpers`. In the example below you can find how to evaluate the `MatrixNet`
that is needed by the HLT1TwoTrackMVA line on all the two particle combinations of a 4 body decay.

```
def getVars(i, j):
 
  combination = str(i)+"_"+str(j)

  Vars = {
  'chi2'   : "RELINFO('"+ relinfo_output +"', 'VERTEX_CHI2_COMB_"+combination+"',0)",
  'fdchi2' : "RELINFO('"+ relinfo_output +"', 'VDCHI2_OWNPV_COMB_"+combination+"',0)",
  'sumpt'  : "RELINFO('"+ relinfo_output +"', 'SUMPT_COMB_"+combination+"',0)",
  'nlt16'  : "RELINFO('"+ relinfo_output +"', 'NLT_OWNPV_COMB_"+combination+"',0)"
  }

  return Vars

from MVADictHelpers import addMatrixnetclassifierTuple
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,2), "HLt1TwoTrackMVAEmulations_1_2", True)
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,3), "HLt1TwoTrackMVAEmulations_1_3", True)
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,4), "HLt1TwoTrackMVAEmulations_1_4", True)
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(2,3), "HLt1TwoTrackMVAEmulations_2_3", True)
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(2,4), "HLt1TwoTrackMVAEmulations_2_4", True)
addMatrixnetclassifierTuple(tuple, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(3,4), "HLt1TwoTrackMVAEmulations_3_4", True)

```

## Additional Documentation

Talks in which the usage of this tool has been shown:

- [Analysis and Software week, May 8th 2019](https://indico.cern.ch/event/811274/contributions/3414033/attachments/1840719/3017872/TriggerEmu_in_TO_Simone.pdf)
- [R1R2 performance WG meeting April 3rd 2019](https://indico.cern.ch/event/810314/contributions/3380333/)

