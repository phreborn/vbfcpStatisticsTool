echo "--- setting up the env ---"
lsetup "lcgenv -p LCG_88 x86_64-slc6-gcc62-opt ROOT";
export CC=/cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/6.2.0/x86_64-slc6/bin/gcc; 
export CXX=/cvmfs/sft.cern.ch/lcg/releases/LCG_88/gcc/6.2.0/x86_64-slc6/bin/g++;

echo "--- setting up eigen ---"
lsetup "lcgenv -p LCG_88 x86_64-slc6-gcc62-opt eigen";


