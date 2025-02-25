#PBS -S /bin/tcsh
#PBS -N TESTING123
#PBS -m abe
#PBS -l select=10:ncpus=24:mpiprocs=24:model=has
#PBS -l walltime=1:00:00
#PBS -k oe
#PBS -r n
#PBS -q devel
#PBS -W group_list=s1152
source /usr/share/modules/init/csh
module load mpi-sgi/mpt comp-intel/2018.3.222 python3/3.7.0
source /usr/local/lib/global.cshrc

setenv MPI_REQUEST_MAX 512
setenv MPI_SHEPHERD true
setenv MPI_BUFS_PER_PROC 512

setenv WDIR /home4/bburning/retrievals/longWaveMie

setenv PATH ${PATH}:${WDIR}:/u/scicon/tools/bin
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}

#setenv OMP_NUM_THREADS 20

unlimit stacksize

limit coredumpsize 0

set time_start=`date '+%T%t%d_%h_06'`
  
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------



cd ${WDIR}


mpiexec -np 240 python TEST_here.py > /nobackup/bburning/TEST.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

