from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.lammps import mk_nvt_input_fr_moltemplate
from AtomicVirtuaLab.espresso import mk_qe_input_vcrelax
import AtomicVirtuaLab.globalv as g
from ase.io import read, write
from ase.visualize import view
import os

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'

smiles={
    'pmma_4mer':'CCC(C)(CC(C)(CC(C)(CC(C)(C)C(=O)OC)C(=O)OC)C(=O)OC)C(=O)OC'
}

mollist={
    'pmma_4mer':1
}

x_box=20.0
y_box=20.0
z_box=20.0

os.makedirs('./pmma_4mer',exist_ok=True)
os.chdir('./pmma_4mer')

for smi in smiles:
    smiles2xyz(smiles[smi],smi,True)

pmma = read('./pmma_4mer.xyz')

#view(pmma)

del pmma[[atom.index for atom in pmma if atom.index in [0,30,31,32,13,53,54,55]]]

shift = pmma[10].position

pmma.translate(-shift)
pmma.translate([0,12.5,12.5])

pmma.set_cell([8.5,25,25])

#pmma = read('./tmp.cif')
view(pmma)

mk_qe_input_vcrelax(pmma,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,kpts=None,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})


mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cp_rdlt()

for smi in smiles:
    mklt(smiles[smi],smi)

mk_system_lt(mollist,x_box,y_box,z_box)

os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')

symbols = get_chemical_symbols('system.data')

mk_nvt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,12345,False)

os.chdir('../')

