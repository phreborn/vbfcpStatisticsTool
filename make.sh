echo "--- cd build ---"
cd build;
echo "--- cmake ---"
cmake .. -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt -DBOOST_INCLUDEDIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt/include/boost -DBoost_LIBRARY_DIRS=/cvmfs/sft.cern.ch/lcg/releases/LCG_93/Boost/1.66.0/x86_64-centos7-gcc62-opt/lib;
echo "--- make ---"
make VERBOSE=1;
echo "--- go back ---"
cd ..
