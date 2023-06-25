from AtomicVirtuaLab.vasp import mk_cellopt_incar, mk_kpoints, mk_potcar, mk_poscar, mk_dos_incar, mk_band_incar, \
                                 plot_vasp_band, plot_vasp_dos
from AtomicVirtuaLab.io import rd_cif, rd_contcar
import AtomicVirtuaLab.globalv as g
import os

g.vasppot = '/media/sf_nanoVM/myVASP/potpaw'
g.cifs = '/media/sf_nanoVM/myPython/AtomicVirtuaLab/cifs'
g.contcars = '/media/sf_nanoVM/myVASP/CuI/cellopt'

#cell = rd_cif(g.cifs+'/'+'CuI_zincblende.cif',primitive_cell=True)
cell = rd_contcar(g.contcars+'/CONTCAR') 

"""
os.makedirs('./vasp/cellopt',exist_ok=True)
os.chdir('./vasp/cellopt')
mk_poscar(cell)
mk_potcar(cell,g.vasppot)
mk_cellopt_incar(cell)
mk_kpoints(4,4,4)
os.chdir('../..')

os.makedirs('./vasp/dos',exist_ok=True)
os.chdir('./vasp/dos')
mk_poscar(cell)
mk_potcar(cell,g.vasppot)
mk_dos_incar(cell)
mk_kpoints(10,10,10)
os.chdir('../..')
os.makedirs('./vasp/band',exist_ok=True)
os.chdir('./vasp/band')
mk_poscar(cell)
mk_potcar(cell,g.vasppot)
mk_band_incar(cell)
mk_kpoints(10,10,10,band=True)
os.chdir('../..')
"""

#plot_vasp_band('/media/sf_nanoVM/myVASP/CuI/band/vasprun.xml')
plot_vasp_dos('/media/sf_nanoVM/myVASP/CuI/dos/vasprun.xml')

