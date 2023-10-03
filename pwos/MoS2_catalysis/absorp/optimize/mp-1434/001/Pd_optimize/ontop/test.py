from ase.io import read,write
from ase.visualize import view

cell = read('./qe_relax.pwo',index=':')

view(cell)

for atoms in cell:
  #atoms.sort()
  atoms.write('test.extxyz',format='extxyz',append=True,)

