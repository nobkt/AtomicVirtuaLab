from AtomicVirtuaLab.io import smiles2xyz
from ase.io import read
from ase.visualize import view
import os


os.makedirs('xyzfiles',exist_ok=True)
os.chdir('xyzfiles')

molname='HSO4-'
smiles = 'OS([O-])(=O)=O'
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)

