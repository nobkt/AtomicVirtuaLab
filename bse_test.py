from AtomicVirtuaLab.espresso import mk_qe2yambo_input_scf, mk_qe2yambo_input_nscf
from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
import os

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

os.makedirs('bse_test',exist_ok=True)
os.chdir('bse_test')
cell = rd_cif(g.cifs+'/'+'PY138_crystal.cif',primitive_cell=False)

os.makedirs('PWSCF',exist_ok=True)
os.chdir('PWSCF')

os.makedirs('input',exist_ok=True)
os.chdir('input')
mk_qe2yambo_input_scf(cell,'pbe','us',level='high',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='auto',tstress=False,nosym=False,options={},nspin=False)
mk_qe2yambo_input_nscf(cell,'pbe','us',level='high',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='auto',tstress=False,nosym=False,options={},nspin=False)



