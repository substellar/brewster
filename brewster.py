#PBS -S /bin/tcsh
#PBS -N D1425_nc
#PBS -m abe
#PBS -l select=1:ncpus=1:mpiprocs=1:model=has+5:ncpus=17:mpiprocs=17:model=has
#PBS -l walltime=00:15:00
#PBS -k oe
#PBS -r n
#PBS -q devel
#PBS -W group_list=s1152
source /usr/share/modules/init/csh
module load mpi-sgi/mpt.2.15r20 comp-intel/2016.2.181 comp-intel/2016.2.181
source /usr/local/lib/global.cshrc

setenv WDIR /home1/bburning/retrievals/ABDor

setenv PATH ${PATH}:${WDIR}:/u/scicon/tools/bin
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}
setenv MPI_SHEPHERD true
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



cd ${WDIR}

#mpdboot --file=$PBS_NODEFILE --ncpus=1 --totalnum=`cat $PBS_NODEFILE | sort -u | wc -l` --ifhn=`head -1 $PBS_NODEFILE` --rsh=ssh --mpd=`which mpd` --ordered

mpiexec -np 86 python D1425_nc.py > /nobackup/bburning/brew_D1425_nc.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

# terminate the MPD daemon

#mpdallexit
