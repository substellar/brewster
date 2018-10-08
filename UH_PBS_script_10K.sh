#PBS -S /bin/tcsh
#PBS -N 2M0559_NC_10K
#PBS -m abe
#PBS -l nodes=5:ppn=32
#PBS -l walltime=01:30:00
#PBS -k oe
#PBS -q main
source ~/.tcshrc
module load use.own
module unload mpich2-x86_64
#module load intel-mpi
module load mpich2-intel


setenv WDIR /home/bb/retrievals/PY3TEST/brewster/

setenv PATH ${PATH}:${WDIR}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}

py3up

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

#mpiexec -machinefile $PBS_NODEFILE -np 116 python 2M0559_NC_10K.py > /beegfs/car/bb/2M0559_NC_10K.log

mpirun -env I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=1 \
       -machinefile $PBS_NODEFILE -n 160 -ppn 32 \
       python 2M0559_NC_10K.py > /beegfs/car/bb/2M0559_NC_10K.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

