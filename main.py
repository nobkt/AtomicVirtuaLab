from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.lammps import mk_npt_compress_input_fr_moltemplate, mk_npt_input_fr_moltemplate
from ase.io import read
import shutil
import os
import sys


os.makedirs('./testbox',exist_ok=True)
os.chdir('./testbox')

smiles={
    'DGEBA':'CC(C)(C1=CC=C(OCC2CO2)C=C1)C1=CC=C(OCC2CO2)C=C1',
    'DDS':'NC1=CC(=CC=C1)S(=O)(=O)C1=CC(N)=CC=C1'
}

mollist={
    'DGEBA':26,
    'DDS':26
}


x_box=50.0
y_box=50.0
z_box=50.0

os.makedirs('./packmol_files',exist_ok=True)
os.chdir('./packmol_files')
packmolpath = os.getcwd()

for smi in smiles:
    smiles2xyz(smiles[smi],smi,True)

mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
os.chdir('../')

os.makedirs('./rdlt_files',exist_ok=True)
os.chdir('./rdlt_files')
rdltpath = os.getcwd()

cp_rdlt()

for smi in smiles:
    mklt(smiles[smi],smi)
os.chdir('../')

os.makedirs('moltemplate_files',exist_ok=True)
os.chdir('moltemplate_files')
moltemplatepath = os.getcwd()
shutil.copy(packmolpath+'/system.xyz','./')
for smi in smiles:
    shutil.copy(rdltpath+'/'+str(smi)+'.lt','./')
mk_system_lt(mollist,x_box,y_box,z_box)
os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')
symbols = get_chemical_symbols('system.data')
os.chdir('../')

os.makedirs('npt_compress',exist_ok=True)
os.chdir('npt_compress')
shutil.copy(moltemplatepath+'/system.data','./')
shutil.copy(moltemplatepath+'/system.in.charges','./')
shutil.copy(moltemplatepath+'/system.in.init','./')
shutil.copy(moltemplatepath+'/system.in.settings','./')
mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,2000,2000,10,10000,400,100,10000,400,2000000,300,1.0,'iso',12345,False,qeq=True)
os.chdir('..')

os.makedirs('npt',exist_ok=True)
os.chdir('npt')
shutil.copy(moltemplatepath+'/system.data','./')
shutil.copy(moltemplatepath+'/system.in.charges','./')
shutil.copy(moltemplatepath+'/system.in.init','./')
shutil.copy(moltemplatepath+'/system.in.settings','./')
mk_npt_input_fr_moltemplate(symbols,False,0.5,1,1,10,300,1.0,'iso',12345,False,qeq=True)
os.chdir('..')

#mk_nvt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,12345,False)
#mk_npt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,0,'iso',12345,False)
#mk_nvtdeform_input_fr_moltemplate(symbols,False,1.0,200,200,20000,300,1.0e10*1.0e-15,200,12345,True)

os.chdir('../')

