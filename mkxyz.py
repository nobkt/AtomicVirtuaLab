from AtomicVirtuaLab.io import smiles2xyz
from ase.io import read
from ase.visualize import view
import os


os.makedirs('xyzfiles',exist_ok=True)
os.chdir('xyzfiles')

molname='isoprene-30'
smiles = 'C\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CC\C=C(/C)CCC=C(C)C'
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)

