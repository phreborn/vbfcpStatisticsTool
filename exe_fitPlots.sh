#!/bin/bash

ds=$(ls ../xmlAnaWSBuilder/run/workspace/ | cut -d _ -f 3 | cut -d . -f 1 | grep -v d)
echo $ds

#for d in ${ds};do
for d in m00;do
  for cat in TT TL LT LL;do
    for bin in b1 b2 b3 b4 b5 b6;do
      bin/readWS.exe ../xmlAnaWSBuilder/run/workspace/vbf_cp_${d}/vbf_cp_${d}.root ${cat}_${bin} ${d} 
    done
  done
done
