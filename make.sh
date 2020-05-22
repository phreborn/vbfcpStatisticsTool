echo "--- cd build ---"
cd build;
echo "--- cmake ---"
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/include/boost -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_97_ATLAS_1/Boost/1.72.0/x86_64-centos7-gcc8-opt/lib
echo "--- make ---"
#make VERBOSE=1;
make
echo "--- go back ---"
cd ..
