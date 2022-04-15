#!/bin/bash

wsDir=wsToy1
wsDir=WSAllCats_SMEFT
wsDir=WSAllCats

ds=$(ls ../xmlAnaWSBuilder/run/${wsDir}/ | cut -d _ -f 3 | cut -d . -f 1 | grep -v d)
echo $ds

dataset=asimovData_SB_SM
dataset=combData

#for d in ${ds};do
#for d in p0d25;do
for d in p01;do
  for cat in TT TL LT;do
    for bin in b1 b2 b3 b4 b5 b6;do
      bin/readWS_fit.exe ../xmlAnaWSBuilder/run/${wsDir}/vbf_cp_${d}/vbf_cp_${d}.root ${cat}_${bin} ${d} ${dataset}
    done
  done
done
