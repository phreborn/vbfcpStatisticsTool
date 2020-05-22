# lxplus environment

Setup the ATLAS environment:

~~~~
setupATLAS;
~~~~

Get packages:

~~~~
lsetup git
git clone https://:@gitlab.cern.ch:8443/mstankai/StatisticsTools.git
cd StatisticsTools
git clone https://:@gitlab.cern.ch:8443/atlas-publications-committee/atlasrootstyle.git
git clone https://github.com/haarcuba/cpp-text-table.git
~~~~

Setup `ROOT` via an `LCG` release:

~~~~
lsetup "views LCG_97_ATLAS_1 x86_64-centos7-gcc8-opt"
export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0/x86_64-centos7/bin/gcc
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/8.3.0/x86_64-centos7/bin/g++
~~~~


# Compiling tools

~~~~
mkdir build;
cd build;
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-op-DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/include/boost -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/include/boost -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/lib

make VERBOSE=1;
~~~~

# Running tools

## Quick fit

~~~~
./bin/run_fit.exe --input input.root --workspace combined --data asimovData --poi mu
~~~~

## Nuisance impact (correlation)

~~~~
./bin/run_pulls.exe --input input.root --workspace combined --data asimovData --poi mu --parameter alpha_sys
~~~~

# Plotting output

## Likelihood scan

### 1D

~~~
./bin/plot_scan.exe --input "root-files/scan/scan" --poi mu --x 0.0 2.0 --y 0.0 20.0 --color "#0000ff" --style 1 --legend myScan --axis_label mu
~~~

### 2D

~~~
./bin/plot_scan.exe --input "root-files/scan/scan" --poi mu1 mu2 --x 0.0 2.0 --y -5.0 5.0 --color "#0000ff" --style 1 --legend myScan --axis_label mu1 mu2
~~~
# Submitting to lxbatch (condor)

## Likelihood scan

### 1D

~~~~
python submit_scan.py input.root --workspaceName combined --data asimovData --poi mu1 --folder scan --scanRange 0.0:2.0 --bins 10 --pointsPerJob 4 --condor
~~~~

### 2D

~~~~
python submit_scan.py input.root --workspaceName combined --data asimovData --poi "mu1,mu2" --folder scan --scanRange "0.0:2.0,-5.0:5.0" --bins "10,10" --pointsPerJob 8 --condor
~~~~

## Nuisance parameter impact (correlation between POI and NP)

~~~~
python submit_pulls.py input.root --workspaceName combined --data asimovData --poi mu --folder ranking --parametersPerJob 2 --condor
~~~~
