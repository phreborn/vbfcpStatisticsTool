# Finds the required directories to include Eigen. Since Eigen is
# only header files, there is no library to locate, and therefore
# no *_LIBRARIES variable is set.

include(FindPackageHandleStandardArgs)
unset(EIGEN_FOUND)

find_path(Eigen_INCLUDE_DIR
        NAMES
        eigen3
        eigen3/unsupported
        eigen3/Eigen
        PATH_SUFFIXES include
        PATHS
        /cvmfs/sft.cern.ch/lcg/releases/LCG_88/eigen/3.2.9/x86_64-slc6-gcc62-opt)

# set Eigen_FOUND
find_package_handle_standard_args(Eigen DEFAULT_MSG Eigen_INCLUDE_DIR)

# set external variables for usage in CMakeLists.txt
if(EIGEN_FOUND)
    set(Eigen_INCLUDE_DIRS ${Eigen_INCLUDE_DIR} ${Eigen_INCLUDE_DIR}/eigen3)
endif()

# hide locals from GUI
mark_as_advanced(Eigen_INCLUDE_DIR)
