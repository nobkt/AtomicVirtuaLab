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

aprun -n 36 -N 36 -j 1 /work/app/QuantumESPRESSO/qe-7.2/bin/pw.x < qe_scf.pwi 1> qe_scf.pwo 2> qe_scf.err
aprun -n 36 -N 36 -j 1 /work/app/QuantumESPRESSO/qe-7.2/bin/dos.x < qe_dos.pwi 1> qe_dos.pwo 2> qe_dos.err
aprun -n 36 -N 36 -j 1 /work/app/QuantumESPRESSO/qe-7.2/bin/projwfc.x < qe_projwfc.pwi 1> qe_projwfc.pwo 2> qe_projwfc.err

cd; if cp -raf $WORKDIR/$DIRNAME $PBS_O_WORKDIR/.. ; then rm -rf $WORKDIR; fi

