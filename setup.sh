#!/bin/bash

echo "--- setting up the env ---"
lsetup "lcgenv -p LCG_93 x86_64-centos7-gcc62-opt ROOT";

echo "--- setting up eigen ---"
lsetup "lcgenv -p LCG_93 x86_64-centos7-gcc62-opt eigen";

echo "--- setting up paths ---"
export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/gcc;
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/g++;

echo "--- setting up cmake ---"
lsetup cmake
