from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.build import make_supercell
from AtomicVirtuaLab.build import slabgen

g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

cell = rd_cif(g.cifdir+'/MoS2_mp1434.cif')

nx=1
ny=1
nz=1

#cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)
cell = slabgen(cell,0,0,1,nx,ny,nz,10.0,10.0)

view(cell)
