#!/bin/bash

AsimovType=floatMu
dataset=asimovData_SM_floatMu
AsimovType=floatMuFixMu1
dataset=asimovData_SM_Mu1

AsimovType=floatMu_fixsyst
dataset=asimovData_SM_floatMu
AsimovType=floatMuFixMu1_fixsyst
dataset=asimovData_SM_Mu1

AsimovType=blind
dataset=asimovData_SB_SM

suffix=_${AsimovType}

for dval in m00
do
  ./bin/combMyyPlots.exe /publicfs/atlas/atlasnew/higgs/hgg/chenhr/vbfcp/WSconfigs/Expected/blind_vs_unblind/${AsimovType}/outAllCats_allSys/out_${dval}.root ${dataset} ${suffix}
done
