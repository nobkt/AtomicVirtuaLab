from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.build import make_supercell
from AtomicVirtuaLab.build import slabgen
from ase.io import read
from AtomicVirtuaLab.lammps import mk_nvt_input_deepmd_friction
import shutil
import os
import sys


g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"


# make cell
g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

cell1 = rd_cif(g.cifdir+'/MoS2_mp1434.cif')
cell2 = rd_cif(g.cifdir+'/MoS2_mp2815.cif')


# 6-layers
#nx=9
#ny=5
#nz1=2
#nz2=3
#nlayer1=6
#nlayer2=6

# 6-layers xy expand
#nx=18
#ny=10
#nz1=2
#nz2=3
#nlayer1=6
#nlayer2=6

# 20-layers
nx=9
ny=5
nz1=7
nz2=10
nlayer1=21
nlayer2=20

os.makedirs('./MoS2_friction',exist_ok=True)
os.chdir('./MoS2_friction')

cell1 = slabgen(cell1,0,0,1,1,1,nz1,10.0,10.0)
os.system('atomsk slab_3.cif -orthogonal-cell MoS2_mp1434_slab_3_ortho.cfg')
#view(cell1)
cell1 = read('MoS2_mp1434_slab_3_ortho.cfg')

cell2 = slabgen(cell2,0,0,1,1,1,nz2,10.0,10.0)
#view(cell2)
os.system('atomsk slab_1.cif -orthogonal-cell MoS2_mp2815_slab_1_ortho.cfg')
cell2 = read('MoS2_mp2815_slab_1_ortho.cfg')

#cell1 = make_supercell(cell1,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
#cell2 = make_supercell(cell2,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
#view(cell1)
#view(cell2)

os.makedirs('./mp1434',exist_ok=True)
os.chdir('./mp1434')
#shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
cell1 = make_supercell(cell1,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
view(cell1)
cell1.write('cell1.cif')
os.makedirs('./y-friction',exist_ok=True)
os.chdir('./y-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell1,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer1,vdirec='y')
#mk_nvt_input_deepmd_friction(cell1,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='y')
os.chdir('../')
os.makedirs('./x-friction',exist_ok=True)
os.chdir('./x-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell1,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer1,vdirec='x')
#mk_nvt_input_deepmd_friction(cell1,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='x')
os.chdir('../')
os.chdir('../')

cell2 = read('./MoS2_mp2815_slab_1_ortho.cfg')

os.makedirs('./mp2815',exist_ok=True)
os.chdir('./mp2815')
cell2 = make_supercell(cell2,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
view(cell2)
cell2.write('cell2.cif')
os.makedirs('./y-friction',exist_ok=True)
os.chdir('./y-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell2,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer2,vdirec='y')
#mk_nvt_input_deepmd_friction(cell2,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='y')
os.chdir('../')
os.makedirs('./x-friction',exist_ok=True)
os.chdir('./x-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell2,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer2,vdirec='x')
#mk_nvt_input_deepmd_friction(cell2,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='x')
os.chdir('../')

print(len(cell1))
print(len(cell2))