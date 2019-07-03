#!/bin/bash

lsetup "lcgenv -p LCG_93 x86_64-centos7-gcc62-opt ROOT";

export CC=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/gcc;
export CXX=/cvmfs/sft.cern.ch/lcg/releases/gcc/6.2.0-b9934/x86_64-centos7/bin/g++;
