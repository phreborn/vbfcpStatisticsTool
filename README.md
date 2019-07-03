# lxplus environment

Setup the ATLAS environment:

~~~~
setupATLAS;
~~~~

Setup `ROOT`, `eigen` and `yaml-cpp` via `LCG` releases:

~~~~
lsetup "lcgenv -p LCG_93 x86_64-centos7-gcc62-opt ROOT";
export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/gcc;
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/g++;
~~~~

Setup a recent `cmake` version:

~~~~
lsetup cmake;
~~~~

# Compiling tools

~~~~
cd StatisticsTools;
mkdir build;
cd build;
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt/include/boost -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt/lib;
make VERBOSE=1;
~~~~

# Running tools

## Quick fit

~~~~
./bin/fit.exe --input input.root --workspace combined --data asimovData --poi mu
~~~~

## Nuisance impact (correlation)

~~~~
./bin/pulls.exe --input input.root --workspace combined --data asimovData --poi mu --parameter alpha_sys
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
