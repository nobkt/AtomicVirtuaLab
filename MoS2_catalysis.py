import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.io import rd_cif, smiles2xyz
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.packmol import mk_packmol_slab_random
from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.siesta import mk_siesta_input_npt
from ase.io import read
from ase.build import make_supercell
from ase.visualize import view
import os
import sys

g.qepot = '/media/sf_nanoVM/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/media/sf_nanoVM/myPython/AtomicVirtuaLab/siesta_pseudo'

g.cifdir='./cifs'
cell = rd_cif(g.cifdir+'/MoS2_mp2815.cif')

os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

# <001> surface
os.makedirs('./001',exist_ok=True)
os.chdir('./001')
slab001 = slabgen(cell,0,0,1,1,1,1,10.0,10.0)
os.system('atomsk slab_1.cif -orthogonal-cell slab_1_ortho.cfg')
slab001 = read('slab_1_ortho.cfg')
slab001 = make_supercell(slab001,([5,0,0],[0,3,0],[0,0,1]),wrap=True)
slab001.write('slab.xyz')
lat = slab001.get_cell()
z=[]
for atom in slab001:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)
#view(slab001)
smiles2xyz('Cl[Pd]Cl','PdCl2',False,userandom=True)
smiles2xyz('O','H2O',True)
smiles2xyz('CCO','CH3CH2OH',True)
mol_PdCl2 = read('PdCl2.xyz')
mol_H2O = read('H2O.xyz')
mol_CH3CH2OH = read('CH3CH2OH.xyz')
#view(mol_PdCl2)
#view(mol_H2O)
#view(mol_CH3CH2OH)
buffer=2.0
slab_top=zmax+buffer
slab_bottom=zmin-buffer
molboxlist={
    'mol1':
        {
            'mol':'PdCl2',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol2':
        {
            'mol':'H2O',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol3':
        {
            'mol':'CH3CH2OH',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol4':
        {
            'mol':'PdCl2',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol5':
        {
            'mol':'H2O',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol6':
        {
            'mol':'CH3CH2OH',
            'num':4,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        }
}
mk_packmol_slab_random('slab',molboxlist,shift=0.5)
os.system('packmol < packmol.inp')

cell = read('system.xyz')
cell.set_cell(lat)

#view(cell)
print(len(cell))
mk_qe_input_npt(cell,'pbe','paw',400,100,0,dt=0.5,level='high',estep=1000,nstep=200000,kpts=None,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
mk_siesta_input_npt(cell,'PBE','SZ',25.0,[1,1,1],473.0,0.0,20000,pseudo_path=g.siesta_pot,MixingWeight=0.1,MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
