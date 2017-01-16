#PBS -S /bin/tcsh
#PBS -N 2M2224_nc
#PBS -m abe
#PBS -l select=1:ncpus=1:mpiprocs=1:model=has+5:ncpus=17:mpiprocs=17:model=has
#PBS -l walltime=35:00:00
#PBS -k oe
#PBS -r n
#PBS -q long
#PBS -W group_list=s1152
source /usr/share/modules/init/csh
module load mpi-intel/4.1.1.036 comp-intel/2015.0.090 python/2.7.10
source /usr/local/lib/global.cshrc


setenv PATH ${PATH}:/home1/bburning/retrievals/2M2224spex:/u/scicon/tools/bin
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home1/bburning/retrievals/2M2224spex
#setenv MPI_BUFS_PER_PROC 512
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



cd /home1/bburning/retrievals/2M2224spex

mpdboot --file=$PBS_NODEFILE --ncpus=1 --totalnum=`cat $PBS_NODEFILE | sort -u | wc -l` --ifhn=`head -1 $PBS_NODEFILE` --rsh=ssh --mpd=`which mpd` --ordered

mpiexec -machinefile $PBS_NODEFILE -np 86 python brewster_nc.py > /nobackup/bburning/brew_2M2224spex_nc.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

# terminate the MPD daemon

mpdallexit
