from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.lammps import mk_nvt_input_uff
from AtomicVirtuaLab.io import rd_cif, rd_lammpsdata
from AtomicVirtuaLab.build import sortmol
from ase.build import make_supercell
from ase.visualize import view
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.tools import add_displacement, deform_cell, scale_cell
from AtomicVirtuaLab.siesta import mk_siesta_input_npt
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'

# FAPbBr3.cif
os.makedirs('./FAPbBr3/data/bulk',exist_ok=True)
os.chdir('./FAPbBr3/data/bulk')

nx=2
ny=2
nz=2

cell = rd_cif(g.cifs+'/'+'FAPbBr3.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)
#print(len(cell))
#sys.exit()
P = [0.0]
T=[100,300,500,700,900,1100]

for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        #cell_liquid = cell.copy()
        #cell_liquid = scale_cell(cell_liquid,0.8)
        #cell_liquid = add_displacement(cell_liquid,0.05)
        #cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #cell_liquid = sortmol(cell_liquid)
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        #mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=0.5,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        mk_siesta_input_npt(cell_solid,'PBEsol','DZP',100.0,None,T0,P0,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
        os.chdir('../')
        #
        # liquid model
        #
        #os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        #os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        #os.makedirs('./tmp',exist_ok=True)
        #os.chdir('./tmp')
        #mk_nvt_input_uff(cell_liquid,0.0001,200,200,20000,100.0,12345)
        #os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        #cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        #view(cell_liquid)
        #sys.exit()
        #os.chdir('../')
        #os.system('rm -rf tmp')
        #mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=0.5,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        #os.chdir('../')


    