#PBS -S /bin/tcsh
#PBS -N YOUR_Runname
#PBS -m abe
#PBS -l nodes=8:ppn=32
#PBS -l walltime=100:00:00
#PBS -k oe
#PBS -q main

source ~/.tcshrc
module unload mpich2-x86_64
module load mpich2-intel

# This should be your working directory
setenv WDIR /home/bb/retrievals/longWaveMie

# Add the working directory to your path
setenv PATH ${PATH}:${WDIR}
# Add WD to LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}

# Active python 3 environment
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


mpirun -env I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=1 \
       -machinefile $PBS_NODEFILE -n 256 -ppn 32 \
       python YOUR_BREWSTER.py > /your/path/for/saving/log



set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

