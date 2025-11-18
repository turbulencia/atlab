# BINS_PATH contains the executables to be tested; needs to be passed externally
#
PRJS=\
     Case01 \
     Case02 \
     Case03 \
     Case06 \
     Case07 \
     Case08 \
     Case09 \
     Case11 \
     Case21 \
     Case25 \
     Case26 \
     Case27 \
     Case31 \
     Case32 \
     Case41 \
     Case51 \
     Case52 \
     Case61 \
     Case62 \
     Case63 \
     Case64 \
     Case65 \
     Case71 \
     Case72 \
     Case73 
#     Case28 \

FILES=avg* grid* flow* scal* tlab.log tlab.err tlab.war tlab.ini.bak dns.out part* traj* main.* field.* figure1.pdf gmon.out

PROCS=4
MPIRUN="mpirun -np  $(PROCS)"

##########################################################################
#   Directives
##########################################################################
SHELL := /bin/bash

check:
	@echo -e "Running selected cases in examples directory."
	@echo -e "Reference incompressible values obtained with FFTW/ESTIMATE mode."
	@(for prj in $(PRJS); do (echo -n "[`date +%T`] Checking $$prj... "; cd $$prj; ../check_run.sh ../../$(BINS_PATH) ) done )

check-mpi:
	@echo -e "Running selected cases in examples directory in parallel mode."
	@echo -e "Reference incompressible values obtained with FFTW/ESTIMATE mode."
	@(for prj in $(PRJS); do (echo -n "[`date +%T`] Checking $$prj... "; cd $$prj; ../check_run.sh ../../$(BINS_PATH) $(MPIRUN)) done )

check-mixed:
	@echo -e "Running selected cases in examples directory in parallel mode with serial preprocessing."
	@echo -e "Reference incompressible values obtained with FFTW/ESTIMATE mode."
	@(for prj in $(PRJS); do (echo -n "[`date +%T`] Checking $$prj... "; cd $$prj; ../check_run.sh ../../$(BINS_PATH) $(MPIRUN) ../../$(BINS_PATH_PP)) done )

##########################################################################
clean: 
	@(for prj in $(PRJS); do (cd $$prj; rm -f $(FILES) ) done )

proper:
	-/bin/rm -rf *~
	@(for prj in $(PRJS); do (cd $$prj; rm -f $(FILES) ) done )

