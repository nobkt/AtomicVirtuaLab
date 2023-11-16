import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.siesta import mk_siesta_input_scf
from ase.io import read, write
from ase.visualize import view
from ase.build import make_supercell, graphene
from ase.cluster.cubic import FaceCenteredCubic
import os

g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

os.makedirs('test_siesta',exist_ok=True)
os.chdir('test_siesta')

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [4, 3, -1]
trunc = FaceCenteredCubic('Pt', surfaces, layers)
trunc.rotate(45, 'x', rotate_cell=True)
trunc.rotate(54.736, 'z', rotate_cell=True)
trunc.rotate(90, 'y', rotate_cell=True)
#view(trunc)

graphene = graphene(formula='C2', a=2.46, size=(12, 12, 1), vacuum=15.0)
#view(graphene)

comgraphene = graphene.get_center_of_mass()
comtrunc = trunc.get_center_of_mass()

ztrunc=[]
zgraphene=[]
for atom in trunc:
    ztrunc.append(atom.position[2])
for atom in graphene:
    zgraphene.append(atom.position[2])
zmintrunc=min(ztrunc)
zmaxgraphene=max(zgraphene)
shift=[comgraphene[0]-comtrunc[0],comgraphene[1]-comtrunc[1],zmaxgraphene+2.5-zmintrunc]

trunc.translate(shift)

adsorp = graphene+trunc
#view(adsorp)

#print(len(adsorp))

mk_siesta_input_scf(adsorp,'BLYP','DZP',100.0,[1,1,1],g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=10,spin='non-polarized')
adsorp.write('model.cif')