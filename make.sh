echo "--- cd build ---"
cd build;
echo "--- cmake ---"
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/boost-1_62 -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_88/Boost/1.62.0/x86_64-slc6-gcc62-opt/lib;
echo "--- make ---"
make;
echo "--- go back ---"
cd ..
