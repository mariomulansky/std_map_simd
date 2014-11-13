# set include directories

INCLUDES = -I$(BOOST_ROOT) -I${NT2_INCLUDE}
#-I${NT2_INCLUDE}


# GCC COMPILER
# change this if you want to cross-compile
# ARCH = native
# # ARCH = core-avx-i

# CXX = g++
# CC = g++
# CXXFLAGS = -O3 -ffast-math -mtune=${ARCH} -march=${ARCH} -std=c++0x -DNDEBUG ${INCLUDES}

# INTEL COMPILER
# change this if you want to cross-compile
ARCH = Host
# ARCH = AVX
# ARCH = SSE4.2

CXX = icpc
CC = icpc
CXXFLAGS = -O3 -x${ARCH} -std=c++0x -fno-alias -inline-forceinline -DNDEBUG ${INCLUDES}
