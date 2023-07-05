from AtomicVirtuaLab.io import rd_cif, smiles2xyz
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.build import make_supercell
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt
from ase.io import read, write
from ase.geometry.analysis import Analysis
import os

g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

cell = rd_cif(g.cifdir+'/cellose_crystal.cif')
#view(cell)
os.makedirs('cellose',exist_ok=True)
os.chdir('./cellose')

smiles={
    'test2mer':'OCC1OC(OC2C(O)C(O)C(OC3C(O)C(O)C(OC4C(O)C(O)COC4CO)OC3CO)OC2CO)C(O)C(O)C1O'
}

for smi in smiles:
    smiles2xyz(smiles[smi],smi,True)

cp_rdlt()

for smi in smiles:
    mklt(smiles[smi],smi)

cell = slabgen(cell,1,0,0,1,1,1,10.0,50.0)
#os.system('atomsk slab_3.cif -orthogonal-cell MoS2_mp1434_slab_3_ortho.cfg')
#view(cell)

cell = read('slab_10.cif')
cell = make_supercell(cell,([1,0,0],[0,1,0],[0,0,1]),wrap=True)
view(cell)
cell.write('test.cif')

ana = Analysis(cell)
for b in ana.unique_bonds:
    print(b)

