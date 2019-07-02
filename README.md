# lxplus environment

Setup the ATLAS environment:

~~~~
cd;
setupATLAS;
~~~~

Setup `ROOT`, `eigen` and `yaml-cpp` via `LCG` releases:

~~~~
cd;
lsetup "lcgenv -p LCG_93 x86_64-centos7-gcc62-opt ROOT";
export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/gcc;
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/g++;
~~~~

Setup a recent `cmake` version:

~~~~
cd;
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

# Submitting to lxbatch

## Nuisance impact (correlation)

~~~~
python submit_pulls.py input.root --workspaceName combined --data asimovData --poi mu --folder ranking
~~~~
