import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.io import rd_cif
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.espresso import mk_qe_input_vcrelax, mk_qe_input_relax
from ase.io import read
from ase.visualize import view
from ase.constraints import FixAtoms
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifdir='./cifs'

# スラブの構造最適化

#mpid = 1221498
#mpid = 571589
#mpid = 1552

ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint


mpid = 1552
nx100 = 2
ny100 = 2
nz100 = 3

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/Mo2C_catalysis/bulk/mp'+str(mpid)+'/opt/qe_vc-relax.pwo')

os.makedirs('./Mo2C_catalysis/slab/mp-'+str(mpid)+'/opt',exist_ok=True)
os.chdir('./Mo2C_catalysis/slab/mp-'+str(mpid)+'/opt')

os.makedirs('./101',exist_ok=True)
os.chdir('./101')

#slab = slabgen(cell,1,0,1,1,1,1,7.5,7.5)
#print(slab[0].get_cell())
#print(len(slab[0]),len(slab[1]))
#view(slab)
#sys.exit()

slab = slabgen(cell,1,0,1,nx100,ny100,nz100,7.5,7.5)

#slab_c = slab[0].copy()
slab_Mo = slab[2].copy()

print(len(slab_Mo))

#view(slab_c)
#view(slab_Mo)
#sys.exit()

#os.makedirs('./C_surface',exist_ok=True)
#os.chdir('./C_surface')


#c = FixAtoms(indices=[atom.index for atom in slab_c if atom.position[2] < 12])
#slab_c.set_constraint(c)
#view(slab_c)

#z=[]
#for atom in slab_c:
#    z.append(atom.position[2])
#zmax=max(z)
#zmin=min(z)

#shift = (zmax + zmin)/2.0

#slab_c.translate([0.0,0.0,-shift])

#mk_qe_input_relax(slab_c,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

#view(slab_c)

#os.chdir('../')

os.makedirs('./Mo_surface',exist_ok=True)
os.chdir('./Mo_surface')


zlow=10.0
c = FixAtoms(indices=[atom.index for atom in slab_Mo if atom.position[2] < zlow])
slab_Mo.set_constraint(c)
#view(slab_Mo)

z=[]
for atom in slab_Mo:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab_Mo.translate([0.0,0.0,-shift])

mk_qe_input_relax(slab_Mo,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

view(slab_Mo)

os.chdir('../')

# スラブの構造最適化 終了



"""
# バルクの構造最適化
ecutwfc=77.0
ecutrho=539.0
k0 = 4

#mpid = 1221498
#mpid = 571589
mpid = 1552


cell = rd_cif(g.cifdir+'/Mo2C_mp-'+str(mpid)+'.cif')
view(cell)

os.makedirs('./Mo2C_catalysis',exist_ok=True)
os.chdir('./Mo2C_catalysis')

os.makedirs('./bulk',exist_ok=True)
os.chdir('./bulk')

os.makedirs('./mp'+str(mpid),exist_ok=True)
os.chdir('./mp'+str(mpid))

os.makedirs('./opt',exist_ok=True)
os.chdir('./opt')
mk_qe_input_vcrelax(cell,'pbe','paw',level='high',nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)

#バルクの構造最適化 終了
"""