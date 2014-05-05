## Intel compiler
CC = icc 
CFLAGS = -lgsl -lgslcblas -openmp -O2 -I $(INCLUDE)

## GCC compiler
# CC = g++-4.8  # (OpenMP support for gcc >= 4.2)
# CFLAGS = -lgsl -framework vecLib -fopenmp -O2

FILES = runimis.cpp imis.cpp likelihood.cpp model.cpp states.cpp

## Note: OpenMP not supported on LLVM 5.0 compiler (default in Mac OS 10.8)
## Note: may need to add -I <path to GSL headers> and -L <path to GSL libraries> if not in default search path
## Note: -lgslcblas causes double free error on Mac OS X 10.8, -framework vecLib works

all: ancprev ancprev-no2012 adultprev adultprev-no2012

ancprev:
	$(CC) $(FILES) $(CFLAGS) -D ancprev=1 -D incl2012=1 -o runimis-ancprev

ancprev-no2012:
	$(CC) $(FILES) $(CFLAGS) -D ancprev=1 -D incl2012=0 -o runimis-ancprev-no2012

adultprev:
	$(CC) $(FILES) $(CFLAGS) -D ancprev=0 -D incl2012=1 -o runimis-adultprev

adultprev-no2012:
	$(CC) $(FILES) $(CFLAGS) -D ancprev=0 -D incl2012=0 -o runimis-adultprev-no2012
