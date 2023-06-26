from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.nwchem import mk_nwchem_input_opt
from ase.io import read
from ase.visualize import view
import os

# optimize
molname='PY129'
smiles = 'OC1=CC=CC=C1\\N=C\\C1=C(O)C=CC2=C1C=CC=C2'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY138'
smiles = 'ClC1=C(Cl)C(Cl)=C(Cl)C2=C1C(=O)N(C2=O)C1=CC=CC2=C1NC(C=C2)=C1C(=O)C2=C(C1=O)C(Cl)=C(Cl)C(Cl)=C2Cl'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY139'
smiles = 'O=C1NC(=O)C(=C2NC(C3=C2C=CC=C3)=C2C(=O)NC(=O)NC2=O)C(=O)N1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY150'
smiles = 'OC1=C(\\N=N\\C2=C(O)NC(=O)NC2=O)C(=O)NC(=O)N1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY185'
smiles = 'CNC(=O)C(\\C#N)=C1/NC(C2=C1C=CC=C2)=C1C(=O)NC(=O)NC1=O'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

"""
molname='PY129_Cu'
smiles = '[#8]-1[Cu]~2[#8]-[#6]-3=[#6](-[#6]=[#7]~2-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-1-2)-[#6]-1=[#6](-[#6]=[#6]-[#6]=[#6]-1)-[#6]=[#6]-3'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True,smarts=True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=2)
os.chdir('../../../')

molname='PY150_Ni'
smiles = 'O=[#6]-1-[#7]-[#6]-2-[#8][Ni]~3[#8]-[#6]-4-[#7]-[#6](=O)-[#7]-[#6](=O)-[#6]-4-[#7]=[#7]~3-[#6]-2-[#6](=O)-[#7]-1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True,smarts=True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')
"""