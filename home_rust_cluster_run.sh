#!/bin/bash
#SBATCH --job-name=vasp_rust
#SBATCH --partition=dev_single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=25:00
###SBATCH --mem=4000
#SBATCH --mem-per-cpu=225

temperature=300
start_temperature=2000


echo " "
echo "### Setting up shell environment ..."
echo " "
# source $HOME/venv3.8/bin/activate
unset LANG; export LC_ALL="C"; export MKL_NUM_THREADS=1; export OMP_NUM_THREADS=1
export USER=${USER:=`logname`}
export MOAB_JOBID=${SLURM_JOB_ID:=`date +%s`}
export MOAB_SUBMITDIR=${SLURM_SUBMIT_DIR:=`pwd`}
export MOAB_JOBNAME=${SLURM_JOB_NAME:=`basename "$0"`}
export MOAB_JOBNAME=$(echo "${SLURM_JOB_NAME}" | sed 's/[^a-zA-Z0-9._-]/_/g')
export MOAB_NODECOUNT=${SLURM_JOB_NUM_NODES:=1}
export MOAB_PROCCOUNT=${SLURM_NTASKS:=1}
#ulimit -s 1000000

echo " "
echo "### Printing basic job infos to stdout ..."
echo " "
echo "START_TIME           = `date +'%y-%m-%d %H:%M:%S %s'`"
echo "HOSTNAME             = ${HOSTNAME}"
echo "USER                 = ${USER}"
echo "MOAB_JOBNAME         = ${SLURM_JOB_NAME}"
echo "MOAB_JOBID           = ${SLURM_JOB_ID}"
echo "MOAB_SUBMITDIR       = ${SLURM_SUBMIT_DIR}"
echo "MOAB_NODECOUNT       = ${SLURM_JOB_NUM_NODES}"
echo "MOAB_PROCCOUNT       = ${SLURM_NTASKS}"
echo "SLURM_NODELIST       = ${SLURM_NODELIST}"
echo "PBS_NODEFILE         = ${SLURM_JOB_NODELIST}"
if test -f "${SLURM_JOB_NODELIST}"; then
  echo "PBS_NODEFILE (begin) ---------------------------------"
  cat "${SLURM_JOB_NODELIST}"
  echo "PBS_NODEFILE (end) -----------------------------------"
fi
echo "---------------- ulimit -a -S ----------------"
ulimit -a -S
echo "---------------- ulimit -a -H ----------------"
ulimit -a -H
echo "----------------------------------------------"



echo " "
module load numlib/mkl/2021.4.0
module load compiler/intel/2021.4.0
module load mpi/impi/2021.4.0

touch out.txt
time ~/rust_mc/target/release/mc -g ~/111-pair -a $1 -i 1000000 -b 2000 -t temperature -c ~/input_cluster/bulk.poscar  > $MOAB_SUBMITDIR/out.txt
exit_code=$?
echo " "
echo "### Cleaning up files ... removing unnecessary scratch files ..."
echo " "
echo "END_TIME             = `date +'%y-%m-%d %H:%M:%S %s'`"

echo "${SLURM_JOB_ID} is complete: on `date +'%y.%m.%d %H:%M:%S'` ${SLURM_SUBMIT_DIR}" >> ~/job.log

echo " "
echo "###$number Exiting with exit code '$exit_code' ..."
echo " "
exit $exit_code
 
