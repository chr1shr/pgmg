# 3D parallel geometric multigrid makefile
#
# Author : Chris H. Rycroft (UW-Madison / LBL)
# Email  : chr@alum.mit.edu
# Date   : April 1st 2021

# Load the common configuration file
include ../config.mk

iflags=
lflags=-L.

probs=problem_simple.o test_problem.o
l_objs=buffer.o common.o geometry.o region.o region_tr.o geometry.o multigrid.o
objs=$(probs) $(l_objs)
src=$(patsubst %.o,%.cc,$(objs))
execs=gs_parallel mg_parallel mg_test gmap_test factors iir_test lapack_test gmplus_test guide_test

all:
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(mpicxx) $(iflags) -MM $(src) >Makefile.dep


include Makefile.dep

libpgmg.a: $(l_objs)
	rm -f libpgmg.a
	ar rs libpgmg.a $^

objects: $(l_objs)
	rm -f libpgmg.a
	ar rs libpgmg.a $^

gs_parallel: gs_parallel.cc problem_simple.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lpgmg problem_simple.o $(lp_lflags)

mg_parallel: mg_parallel.cc problem_simple.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lpgmg problem_simple.o $(lp_lflags)

mg_test: mg_test.cc problem_simple.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lpgmg problem_simple.o $(lp_lflags)

factors: factors.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

iir_test: iir_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

gmap_test: gmap_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $<

lapack_test: lapack_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lp_lflags)

gmplus_test: gmplus_test.cc problem_simple.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lpgmg problem_simple.o $(lp_lflags)

guide_test: guide_test.cc test_problem.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) -lpgmg test_problem.o $(lp_lflags)

eigenvalues: eigenvalues.cc problem_simple.o libpgmg.a
	$(mpicxx) $(cflags) $(iflags) -o $@ $< $(lflags) $(eigen_iflags) -lpgmg problem_simple.o $(lp_lflags)

%.o: %.cc
	$(mpicxx) -Winline $(cflags) $(iflags) -c $<

clean:
	rm -rf $(execs) $(objs) libpgmg.a *.dSYM

.PHONY: clean all executables depend objects
