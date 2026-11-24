#!/bin/bash
export LC_ALL=en_US.utf8
export PATH=/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr/bin:/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr_sw/bin:/usr/sw/mpi/mpi_20220608_SEA/bin:$PATH
export OMP_NUM_THREADS=1
# Full model blocked in PIO error-path MPI_Reduce with hardware collectives.
# Job 9172322 passed the same point with only this environment setting changed.
export SWCH_USE_HW_COLL=0
