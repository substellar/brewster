#PBS -S /bin/tcsh
#PBS -N D1425_nc
#PBS -m abe
#PBS -l nodes=16:ppn=6
#PBS -l walltime=00:30:00
#PBS -k oe
#PBS -q car

module unload mpich2-x86_64
module load mpich2-intel
#module load mvapich2

#setenv OMP_NUM_THREADS  12

setenv WDIR /home/bb/retrievals/ABDor
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${WDIR}
setenv PATH ${WDIR}:${PATH}    

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


set NPROCS = `wc -l < $PBS_NODEFILE`

echo PBS: NProcs = $NPROCS

cd ${WDIR}


mpiexec -np $NPROCS python D1425_nc_UH.py > /beegfs/car/bb/brew_D1425_nc.log

set time_end=`date '+%T%t%d_%h_06'`
echo Started at: $time_start
echo Ended at: $time_end
echo ------------------------------------------------------
echo Job ends
