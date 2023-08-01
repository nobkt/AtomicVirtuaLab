from ase.build import molecule
from ase.io import read, write
from ase.visualize import view
from AtomicVirtuaLab.espresso import mk_qe_input_relax
import AtomicVirtuaLab.globalv as g
import os

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'

# H2O test
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
h2o_mol = molecule('H2O')
com = h2o_mol.get_center_of_mass()
shift = [5.0,5.0,5.0] - com
h2o_mol.translate(shift)
h2o_mol.set_cell([10.0,10.0,10.0])
#view(h2o_mol)

mk_qe_input_relax(h2o_mol,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,mixing_beta=0.2,kpts=None,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
