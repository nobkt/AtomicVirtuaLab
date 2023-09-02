#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1
#PBS -q P_016
#PBS -N vasp

ulimit -sunlimited

DIRNAME=`basename $PBS_O_WORKDIR`
WORKDIR=/work/scratch/$USER/$PBS_JOBID
mkdir -p $WORKDIR
cp -raf  $PBS_O_WORKDIR $WORKDIR
cd $WORKDIR/$DIRNAME

aprun -n 36 -N 36 -j 1 /work/app/QuantumESPRESSO/qe-7.2/bin/pp.x < pwscf.pp.chg.in 1> pwscf.pp.chg.log 2> pwscf.pp.chg.err
aprun -n 36 -N 36 -j 1 /work/app/QuantumESPRESSO/qe-7.2/bin/pp.x < pwscf.pp.ldos.in 1> pwscf.pp.ldos.log 2> pwscf.pp.ldos.err

cd; if cp -raf $WORKDIR/$DIRNAME $PBS_O_WORKDIR/.. ; then rm -rf $WORKDIR; fi

