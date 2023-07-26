from AtomicVirtuaLab.io import smiles2molandpdb
from ase.io import read
from ase.visualize import view

molname='test'

smiles = 'CC'

mol = smiles2molandpdb(smiles,molname,True)

