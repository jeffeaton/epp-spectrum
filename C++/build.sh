#!/bin/bash

mkfl="CXX=g++-5.1.0 F77=gfortran-5.1.0 PKG_CPPFLAGS=-fpermissive"  # for desktop
## mkfl="CXX=icc F77=ifort"  # for cluster

rm *.o;
echo $mkfl
MAKEFLAGS=$mkfl R CMD SHLIB -lgsl -lgslcblas -lgfortran rlib.cpp model.cpp states.cpp incidence.cpp mvndstpack.f parameters.cpp
