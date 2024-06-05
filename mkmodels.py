from pymatgen.core.surface import generate_all_slabs
from pymatgen.io.ase import AseAtomsAdaptor
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import slabgen
from ase.io import read
from ase.build import make_supercell
from ase.visualize import view
import os
import sys


os.makedirs('mkmodels',exist_ok=True)
os.chdir('mkmodels')
g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/dpdatas'
cell0 = read(g.cifdir+'/MgAl2O4.cif',primitive_cell=True, subtrans_included=False)
cell1 = read(g.cifdir+'/MgAl2O4.cif')
#view(cell1)
#cell0 = make_supercell(cell0,([2,0,0],[0,2,0],[0,0,2]))
cell0.write(g.cifdir+'/MgAl2O4_primitive.cif')
#print(len(cell0))
#view(cell0)

cell = AseAtomsAdaptor.get_structure(cell1)
all_slabs = generate_all_slabs(cell,1,1,15,center_slab=True,primitive=True,in_unit_planes=True,repair=True)
f = open('miller_index.txt','w')
f.write(str(len(all_slabs))+'\n')
print(len(all_slabs))
for slab in all_slabs:
    print(slab.miller_index)
    f.write(str(slab.miller_index)+'\n')
f.close()

slab1 = slabgen(cell1,1,1,1,2,2,2,7.5,7.5)
#view(slab1)
j = 1
for i in [0,2,4,6]:
    slab1[i].write(g.cifdir+'/MgAl2O4_slab111_'+str(j)+'.cif')
    j = j + 1

slab2 = slabgen(cell1,1,1,0,2,1,2,7.5,7.5)
#view(slab2)
j = 1
for i in [0,2]:
    slab2[i].write(g.cifdir+'/MgAl2O4_slab110_'+str(j)+'.cif')
    j = j + 1


slab3 = slabgen(cell1,1,0,0,2,2,1,7.5,7.5)
#view(slab3)
j = 1
for i in [0,2]:
    slab3[i].write(g.cifdir+'/MgAl2O4_slab100_'+str(j)+'.cif')
    j = j + 1
#sys.exit()


"""
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.build import make_supercell
    from ase.io import write
    from ase.visualize import view
    cell = AseAtomsAdaptor.get_structure(cell0)
    slabs = SlabGenerator(cell,(h,k,l),nz,10,in_unit_planes=True).get_slabs(repair=True)
"""