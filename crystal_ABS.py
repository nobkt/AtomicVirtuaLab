from AtomicVirtuaLab.espresso import mk_qe2yambo_input_scf, mk_qe2yambo_input_scf_gaupbe, mk_qe_input_relax
from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.io import read
from ase.build import make_supercell
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/mol_crystals'

os.makedirs('crystal_ABS',exist_ok=True)
os.chdir('crystal_ABS')

l_mols = {
    'PY':['Pig_yellow138_exp_2056_exp_2057'],
    'BQ':['BQ11_exp_2054_exp_2052_final'],
    'QA':['alpha-1','alpha-2','beta','gamma'],
    'perylene':['1140279-pr179','oa1058_structure_1_of_1-pb31','pr149']
}

# YAMBO SCF
os.makedirs('pwscf',exist_ok=True)
os.chdir('pwscf')

ecutwfc=36.0
kx=2
ky=2
kz=2

for mols in l_mols:
    os.makedirs(mols,exist_ok=True)
    os.chdir(mols)
    for mol_ in l_mols[mols]:
        os.makedirs(mol_,exist_ok=True)
        os.chdir(mol_)
        os.makedirs('scf',exist_ok=True)
        os.chdir('scf')
        cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/crystal_ABS/'+str(mols)+'/'+str(mol_)+'/opt/qe_relax.pwo')
        os.makedirs('pbe',exist_ok=True)
        os.chdir('pbe')
        mk_qe2yambo_input_scf(cell,'pbe','nc',level='SG15',estep=9999,ecutwfc=ecutwfc,ecutrho=ecutwfc*4,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='manual',tstress=False,nosym=True,options={},nspin=False)
        os.chdir('../')
        os.makedirs('gaupbe',exist_ok=True)
        os.chdir('gaupbe')
        mk_qe2yambo_input_scf_gaupbe(cell,'pbe','nc',level='SG15',estep=9999,ecutwfc=ecutwfc,ecutrho=ecutwfc*4,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='manual',tstress=False,nosym=True,options={},nspin=False)
        os.chdir('../')
        os.chdir('../')
        os.chdir('../')
    os.chdir('../')

"""
# Optimize
ecutwfc=36
ecutrho=400
kx=2
ky=2
kz=2

for mols in l_mols:
    os.makedirs(mols,exist_ok=True)
    os.chdir(mols)
    for mol_ in l_mols[mols]:
        os.makedirs(mol_,exist_ok=True)
        os.chdir(mol_)
        os.makedirs('opt',exist_ok=True)
        os.chdir('opt')
        cell = rd_cif(g.cifs+'/'+mol_+'.cif',primitive_cell=False)
        mk_qe_input_relax(cell,'pbe','paw',level='high',estep=9999,nstep=9999,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
        os.chdir('../')
        os.chdir('../')
    os.chdir('../')
"""
