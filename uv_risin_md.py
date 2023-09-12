from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mklt_fr_mol, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_npt_input_fr_moltemplate, mk_npt_compress_input_fr_moltemplate
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.io import read
from rdkit import Chem
from rdkit.Chem import Draw, rdDistGeom
import os
import shutil
import sys

g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

nmol = 64

os.makedirs('./uv_risin_md',exist_ok=True)
os.chdir('./uv_risin_md')

mols = Chem.SDMolSupplier(g.cifs+'/uv_risin_molecules.sdf')
i = 1
for mol in mols:
    if i==7:
        rnd = True
    else:
        rnd = True
    os.makedirs('mol'+str(i),exist_ok=True)
    os.chdir('mol'+str(i))
    smiles=Chem.MolToSmiles(mol)
    print('mol'+str(i)+' : ',smiles)
    cp_rdlt()
    mklt(smiles,'mol'+str(i),random=rnd)
    #smiles2xyz(smiles,'mol'+str(i),True,smarts=False,userandom=False)
    Draw.MolToFile(mol,'mol'+str(i)+'.png')
    mollist={
        'mol'+str(i):nmol
    }
    x_box=200.0
    y_box=200.0
    z_box=200.0
    os.makedirs('./packmol',exist_ok=True)
    os.chdir('./packmol')
    smiles2xyz(smiles,'mol'+str(i),True,smarts=False,userandom=rnd)
    view(read('mol'+str(i)+'.xyz'))
    mk_packmol_random(mollist,x_box,y_box,z_box)
    os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
    #tmp = read('system.xyz')
    #tmp.set_cell([x_box,y_box,z_box])
    #view(tmp)
    #tmp.write('test.cif')
    packmolpath = os.getcwd()
    os.chdir('../')
    os.makedirs('moltemplate',exist_ok=True)
    os.chdir('moltemplate')
    shutil.copy('../packmol/system.xyz','./')
    shutil.copy('../mol'+str(i)+'.lt','./')
    mk_system_lt(mollist,x_box,y_box,z_box)
    os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
    os.system('cleanup_moltemplate.sh')
    symbols = get_chemical_symbols('system.data')
    os.chdir('../')
    os.makedirs('npt_p100',exist_ok=True)
    os.chdir('npt_p100')
    shutil.copy('../moltemplate/system.data','./')
    shutil.copy('../moltemplate/system.in.charges','./')
    shutil.copy('../moltemplate/system.in.init','./')
    shutil.copy('../moltemplate/system.in.settings','./')
    #mk_npt_input_fr_moltemplate(symbols,True,0.5,200,200,200000,400,100,'iso',12345,False,qeq=True)
    mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,200,200,10,10000,400,100,10000,400,2000000,300,1.0,'iso',12345,False,qeq=True)
    os.chdir('../')
    i+=1
    os.chdir('../')


