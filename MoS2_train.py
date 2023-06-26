from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.lammps import mk_nvt_input_uff
from AtomicVirtuaLab.io import rd_cif, rd_lammpsdata
from ase.build import make_supercell
from ase.visualize import view
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.tools import add_displacement, deform_cell, scale_cell
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

# mp2815
os.makedirs('./MoS2/data/mp-2815/bulk',exist_ok=True)
os.chdir('./MoS2/data/mp-2815/bulk')

nx=3
ny=3
nz=2

cell = rd_cif(g.cifs+'/'+'MoS2_mp2815.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)
#view(cell)
#print(len(cell))
#sys.exit()
P0 = 0.0
T=[300,1000,2000,3000]

for T0 in T:
    cell_solid = cell.copy()
    cell_solid = add_displacement(cell_solid,0.05)
    cell_solid = deform_cell(cell_solid,0.05,0.5)
    cell_liquid = cell.copy()
    cell_liquid = scale_cell(cell_liquid,0.8)
    cell_liquid = add_displacement(cell_liquid,0.05)
    cell_liquid = deform_cell(cell_liquid,0.05,0.5)
    #
    # solid model
    #
    os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
    os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
    mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
    os.chdir('../')
    #
    # liquid model
    #
    os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
    os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
    os.makedirs('./tmp',exist_ok=True)
    os.chdir('./tmp')
    mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,200000.0,12345)
    os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
    rd_lammpsdata(cell_liquid,'./result.data',True)
    os.chdir('../')
    os.system('rm -rf tmp')
    mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
    os.chdir('../')

    
    
    
    
    

    
    
    
