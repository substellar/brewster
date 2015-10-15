#PBS -S /bin/tcsh
#PBS -N bb_5g_fulltest
#PBS -m abe
#PBS -l select=2:ncpus=21:mpiprocs=21:model=has
#PBS -l walltime=110:00:00
#PBS -k oe
#PBS -q long
#PBS -W group_list=s1152
source /usr/share/modules/init/csh
module load comp-intel/2015.3.187 mpi-mvapich2/1.4.1/intel python/2.7.10
source /usr/local/lib/global.cshrc


setenv PATH ${PATH}:/home1/bburning/retrievals/mycode:/u/scicon/tools/bin
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home1/bburning/retrievals/mycode
setenv OMP_NUM_THREADS 42
setenv MPI_BUFS_PER_PROC 256



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



cd /home1/bburning/retrievals/mycode


mpiexec python brewster.py > brew_fulltest.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends
