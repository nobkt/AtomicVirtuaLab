from AtomicVirtuaLab.vasp import mk_cellopt_incar, mk_kpoints, mk_potcar, mk_poscar, mk_dos_incar, mk_band_incar, \
                                 plot_vasp_band, plot_vasp_dos
from AtomicVirtuaLab.io import rd_cif, rd_contcar
import AtomicVirtuaLab.globalv as g
import os

g.vasppot = '/home/A23321P/modules/applications/vasp/6.4.0/vasppot'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.contcars = '/media/sf_nanoVM/myVASP/CuI/cellopt'

cell = rd_cif(g.cifs+'/'+'I.cif',primitive_cell=False)
#cell = rd_contcar(g.contcars+'/CONTCAR') 


os.makedirs('./vasp/I/cellopt',exist_ok=True)
os.chdir('./vasp/I/cellopt')
mk_poscar(cell)
mk_potcar(cell,g.vasppot)
mk_cellopt_incar(cell)
mk_kpoints(4,4,4)
os.chdir('../../..')

"""
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
"""
plot_vasp_band('/home/A23321P/work/myVASP/CuI/zincblende/band/vasprun.xml')
plot_vasp_dos('/home/A23321P/work/myVASP/CuI/zincblende/dos/vasprun.xml')
"""
