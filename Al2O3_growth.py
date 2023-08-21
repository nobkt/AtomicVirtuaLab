from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.io import rd_cif
from ase.build import make_supercell
from ase.visualize import view
from ase.units import mol
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.tools import add_displacement, deform_cell, scale_cell
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

# FAPbBr3.cif
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')

os.makedirs('a-Al2O3',exist_ok=True)
os.chdir('a-Al2O3')

ldens = 2.93

nx=2
ny=2
nz=1

cell = rd_cif(g.cifs+'/'+'Al2O3.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)

lat = cell.get_cell_lengths_and_angles()
volume = cell.get_volume()
masses = cell.get_masses()
mass = sum(masses)
density = (mass/mol)/(volume*1.0e-24)
scale = density/ldens
scale = (scale)**(1/3)

lat[0] = lat[0]*scale
lat[1] = lat[1]*scale
lat[2] = lat[2]*scale

P=[0]
T=[100,500,1000,2000,3000,4000,5000]

for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        #view(cell_solid)
        #sys.exit()
        cell_liquid = cell.copy()
        cell_liquid.set_cell(lat,scale_atoms=True)
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #view(cell_liquid)
        #sys.exit()
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')

