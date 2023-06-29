from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.build import make_supercell
from AtomicVirtuaLab.build import slabgen
from ase.io import read
import os

g.cifdir='/media/sf_nanoVM/myPython/AtomicVirtuaLab/cifs'

cell1 = rd_cif(g.cifdir+'/MoS2_mp1434.cif')
cell2 = rd_cif(g.cifdir+'/MoS2_mp2815.cif')

nx=18
ny=10
nz1=2
nz2=3

cell1 = slabgen(cell1,0,0,1,1,1,nz1,10.0,10.0)
os.system('atomsk slab_3.cif -orthogonal-cell MoS2_mp1434_slab_3_ortho.cfg')
#view(cell1)
cell1 = read('MoS2_mp1434_slab_3_ortho.cfg')

cell2 = slabgen(cell2,0,0,1,1,1,nz2,10.0,10.0)
#view(cell2)
os.system('atomsk slab_1.cif -orthogonal-cell MoS2_mp2815_slab_1_ortho.cfg')
cell2 = read('MoS2_mp2815_slab_1_ortho.cfg')

cell1 = make_supercell(cell1,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
cell2 = make_supercell(cell2,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
view(cell1)
view(cell2)

