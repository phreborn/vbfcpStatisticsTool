tag=

allJobs=jobsSub.sh
> ${allJobs}

cats=($(cat ../../nom_WS/cats.cfg | grep -v "#" | grep ":" | cut -d : -f 1))

sysList=($(cat ../xmlAnaWSBuilder/run/config/vbf_cp_m00/channel/category_OO_TT_b1.xml | grep Systematic | cut -d '"' -f 2 | grep -v ':category:' | sort | uniq))

for cat in ${cats[*]};do
  sysList[${#sysList[@]}]=ATLAS_Hgg_BIAS_OO_${cat}
done

for sys in ${sysList[@]};do echo ${sys};done

sequence=($(seq 1 10 ${#sysList[@]}))

intvl=9
for init in ${sequence[@]};do
  init=$((${init} - 1))
  fin=$((${init} + ${intvl}))
  jobName=Collect_${init}_${fin}; echo ${jobName}
  #if [ ! -d csv/${jobName} ];then mkdir -p csv/${jobName};fi
  if [ ! -d hep_sub_${jobName} ]; then mkdir hep_sub_${jobName}; fi
  executable=exe_${jobName}.sh
  > ${executable}

  echo "#!/bin/bash" >> exe_${jobName}.sh
  echo "" >> exe_${jobName}.sh
  echo "cd /scratchfs/atlas/huirun/atlaswork/VBF_CP/WSBuilder/StatisticsTools" >> exe_${jobName}.sh
  echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> exe_${jobName}.sh
  echo "source \${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" >> exe_${jobName}.sh
  echo "lsetup \"views LCG_97_ATLAS_1 x86_64-centos7-gcc8-opt\"" >> exe_${jobName}.sh
  echo "export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0/x86_64-centos7/bin/gcc" >> exe_${jobName}.sh
  echo "export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0/x86_64-centos7/bin/g++" >> exe_${jobName}.sh
  for num in `seq ${init} 1 ${fin}`;do
  echo "" >> exe_${jobName}.sh
    #echo "./bin/run_pulls.exe --input ../xmlAnaWSBuilder/run/WSAllCats_SMEFT/vbf_cp_m0d00/vbf_cp_m0d00.root --workspace combWS --data asimovData_SB_SM --poi mu_VBF_RW[0:5] --fix mu_VBF_SM[0],mu[1],mu_ggH[1],mu_ggH_SM[0] --parameter ${sysList[${num}]}" >> exe_${jobName}.sh
    echo "./bin/run_pulls.exe --input ../xmlAnaWSBuilder/run/WSAllCats/vbf_cp_m00/vbf_cp_m00.root --workspace combWS --data asimovData_SB_SM --poi mu_VBF_RW[0:5] --fix mu_VBF_SM[0],mu[1],mu_ggH[1],mu_ggH_SM[0] --parameter ${sysList[${num}]}" >> exe_${jobName}.sh
  done

  chmod +x exe_${jobName}.sh

  echo "hep_sub exe_${jobName}.sh -g atlas -os CentOS7 -wt mid -mem 2048 -o hep_sub_${jobName}/log-0.out -e hep_sub_${jobName}/log-0.err" >> ${allJobs}

  if [ "$(ls hep_sub_${jobName}/)" != "" ];then rm hep_sub_${jobName}/*;fi
done
