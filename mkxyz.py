from AtomicVirtuaLab.io import smiles2xyz
from ase.io import read
from ase.visualize import view

molname='PY129'

smiles = 'OC1=CC=CC=C1\\N=C\\C1=C(O)C=CC2=C1C=CC=C2'
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)

