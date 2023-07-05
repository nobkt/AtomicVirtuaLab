from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.lammps import mk_nvt_input_fr_moltemplate, mk_npt_input_fr_moltemplate,mk_nvtdeform_input_fr_moltemplate
import os

smiles={
    'test2mer':'OCC1OC(OC2C(O)C(O)COC2CO)C(O)C(O)C1O'
}

mollist={
    'test2mer':1
}

x_box=1000.0
y_box=1000.0
z_box=1000.0



os.makedirs('./packmol',exist_ok=True)
os.chdir('./packmol')

for smi in smiles:
    smiles2xyz(smiles[smi],smi,True)

mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cp_rdlt()

for smi in smiles:
    mklt(smiles[smi],smi)

mk_system_lt(mollist,x_box,y_box,z_box)

os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')

symbols = get_chemical_symbols('system.data')

#mk_nvt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,12345,False)
mk_npt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,0,'iso',12345,False)
mk_nvtdeform_input_fr_moltemplate(symbols,False,1.0,200,200,20000,300,1.0e10*1.0e-15,200,12345,True)

os.chdir('../')

