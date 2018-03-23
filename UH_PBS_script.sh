#PBS -S /bin/tcsh
#PBS -N 2M2224_FeEnst2c2P_multi
#PBS -m abe
#PBS -l nodes=6:ppn=18
#PBS -l walltime=80:00:00
#PBS -k oe
#PBS -q main

module load use.own
module unload mpich2-x86_64
module load intel-mpi


setenv WDIR /home/bb/retrievals/longWaveMie

setenv PATH ${PATH}:${WDIR}
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}

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

mpdboot --file=$PBS_NODEFILE --ncpus=1 --totalnum=`cat $PBS_NODEFILE | sort -u | wc -l` --ifhn=`head -1 $PBS_NODEFILE` --rsh=ssh --mpd=`which mpd` --ordered

mpiexec -machinefile $PBS_NODEFILE -np 108 python 2m2224_FeEnst2c2p_multi.py > /beegfs/car/bb/brew_2M2224_FeEnst2c2p.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends

