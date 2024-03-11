#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1
#PBS -q C_002
#PBS -N vasp

module add qe/7.3
export OMP_NUM_THREADS=1

ulimit -sunlimited

DIRNAME=`basename $PBS_O_WORKDIR`
WORKDIR=/work/scratch/$USER/qe/$PBS_JOBID
mkdir -p $WORKDIR
cp -raf  $PBS_O_WORKDIR $WORKDIR
cd $WORKDIR/$DIRNAME

mpirun -np 36 -ppn 36 -hostfile $PBS_NODEFILE pw.x < qe_relax.pwi 1> qe_relax.pwo 2> qe_relax.err
cp ./outdir/pwscf.esm1 ./
rm -rf outdir

cd; if cp -raf $WORKDIR/$DIRNAME $PBS_O_WORKDIR/.. ; then rm -rf $WORKDIR; fi

