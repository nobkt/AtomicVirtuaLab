from AtomicVirtuaLab.io import smiles2xyz
from ase.io import read
from ase.visualize import view

molname='test'

smiles = 'COC1=CC(OC)=CC(OC)=C1'
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)

