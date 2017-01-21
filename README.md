# lxplus environment

Setup the basic environment as follows:

~~~~
cd;
lsetup "lcgenv -p LCG_86 x86_64-slc6-gcc62-opt ROOT";
export CC=/cvmfs/sft.cern.ch/lcg/releases/LCG_86/gcc/6.2.0/x86_64-slc6/bin/gcc; export CXX=/cvmfs/sft.cern.ch/lcg/releases/LCG_86/gcc/6.2.0/x86_64-slc6/bin/g++;
~~~~

# Compiling tools

~~~~
mkdir build;
cd build;
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_86/Boost/1.62.0/x86_64-slc6-gcc62-opt -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_86/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/boost-1_62 -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_86/Boost/1.62.0/x86_64-slc6-gcc62-opt/lib;
make;
~~~~

# Running tools

## Quick fit

~~~~
./bin/fit.exe --input input.root --workspace combined --data asimovData --poi mu
~~~~
