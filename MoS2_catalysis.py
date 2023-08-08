import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.io import rd_cif, smiles2xyz
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.packmol import mk_packmol_slab_random
from AtomicVirtuaLab.espresso import mk_qe_input_npt, mk_qe_input_relax, mk_qe_input_vcrelax, mk_qe_input_dos, mk_qe_input_band,mk_qe_input_scf,plot_qe_dos
from AtomicVirtuaLab.siesta import mk_siesta_input_npt, mk_siesta_input_scf_withEfield
from ase.io import read, write
from ase.build import make_supercell, add_adsorbate
from ase.visualize import view
from ase import Atom
from ase.constraints import FixAtoms
import os
import sys

os.system('export DISPLAY=localhost:11.0')

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'

g.cifdir='./cifs'

"""
# スラブモデルDOS
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

path = os.getcwd()

# mp-1434
os.chdir(path)
mpid = 1434
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')
view(slab)

os.makedirs('./MoS2_catalysis/slab/dos/mp-'+str(mpid)+'/'+str(mill),exist_ok=True)
os.chdir('./MoS2_catalysis/slab/dos/mp-'+str(mpid)+'/'+str(mill))
mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
mk_qe_input_dos(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

# mp-2815
os.chdir(path)
mpid = 2815
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')
view(slab)

os.makedirs('./MoS2_catalysis/slab/dos/mp-'+str(mpid)+'/'+str(mill),exist_ok=True)
os.chdir('./MoS2_catalysis/slab/dos/mp-'+str(mpid)+'/'+str(mill))
mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
mk_qe_input_dos(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

# スラブモデルDOS終了
"""



# 吸着モデル
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]

elm = 'Pd'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_potential',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_potential')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)

# ontop
rlist = [2.1,2.2,2.3,2.4]
adtype = 'ontop'
x0 = 1.577
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge
rlist = [1.0,1.1,1.2,1.3,1.4]
adtype = 'bridge'
x0 = (1.577+4.731)/2.0
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# fcc
adtype = 'fcc'
x0 = 3.154
y0 = 5.463
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# hcp
adtype = 'hcp'
x0 = 1.577
y0 = 6.373
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_potential',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_potential')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# ontop
adtype = 'ontop'
x0 = 1.576
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge
adtype = 'bridge'
x0 = (1.576+4.727)/2.0
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# fcc
adtype = 'fcc'
x0 = (1.576+4.727+3.151)/3.0
y0 = (4.548+4.548+7.278)/3.0
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# hcp
adtype = 'hcp'
x0 = 1.576
y0 = 6.368
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')
# 吸着モデル終了


"""
# MoS2スラブモデル作成
ecutwfc0=77.0
ecutrho0=539.0
kpoint=4
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint


mpid = 1434
nx001 = 4
ny001 = 4
nz001 = 1
nx100 = 5
ny100 = 1
nz100 = 4

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/bulk/optimize/mp'+str(mpid)+'/qe_vc-relax.pwo')

os.makedirs('./MoS2_catalysis/slab/optimize/mp-'+str(mpid),exist_ok=True)
os.chdir('./MoS2_catalysis/slab/optimize/mp-'+str(mpid))

os.makedirs('./001',exist_ok=True)
os.chdir('./001')
#slab = slabgen(cell,0,0,1,1,1,1,7.5,7.5)
#view(slab)
slab = slabgen(cell,0,0,1,nx001,ny001,nz001,7.5,7.5)
id = 2
#view(slab[id])

c = FixAtoms(indices=[atom.index for atom in slab[id] if atom.position[2] < 12])
slab[id].set_constraint(c)

#view(slab[id])

z=[]
for atom in slab[id]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[id].translate([0.0,0.0,-shift])

mk_qe_input_relax(slab[id],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

view(slab[id])

os.chdir('../')


os.makedirs('./100',exist_ok=True)
os.chdir('./100')
slab = slabgen(cell,1,0,0,1,1,1,7.5,7.5)
#view(slab)
#sys.exit()
slab = slabgen(cell,1,0,0,nx100,ny100,nz100,7.5,7.5)
id = 0
#view(slab[0])
#view(slab[1])

c = FixAtoms(indices=[atom.index for atom in slab[0] if atom.position[2] < 9.5])
slab[0].set_constraint(c)

c = FixAtoms(indices=[atom.index for atom in slab[1] if atom.position[2] < 9.5])
slab[1].set_constraint(c)

#view(slab[0])
#view(slab[1])

z=[]
for atom in slab[0]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[0].translate([0.0,0.0,-shift])

z=[]
for atom in slab[1]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[1].translate([0.0,0.0,-shift])

os.makedirs('surface0',exist_ok=True)
os.chdir('surface0')
mk_qe_input_relax(slab[0],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')
os.makedirs('surface1',exist_ok=True)
os.chdir('surface1')
mk_qe_input_relax(slab[1],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

view(slab[0])
view(slab[1])

os.chdir('../')



mpid = 2815
nx001 = 4
ny001 = 4
nz001 = 2
nx100 = 5
ny100 = 2
nz100 = 4

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/bulk/optimize/mp'+str(mpid)+'/qe_vc-relax.pwo')

os.makedirs('./MoS2_catalysis/slab/optimize/mp-'+str(mpid),exist_ok=True)
os.chdir('./MoS2_catalysis/slab/optimize/mp-'+str(mpid))

os.makedirs('./001',exist_ok=True)
os.chdir('./001')
#slab = slabgen(cell,0,0,1,1,1,1,7.5,7.5)
#view(slab)
slab = slabgen(cell,0,0,1,nx001,ny001,nz001,7.5,7.5)
id = 2
#view(slab)
#view(slab[id])

c = FixAtoms(indices=[atom.index for atom in slab[id] if atom.position[2] < 12])
slab[id].set_constraint(c)

#view(slab[id])

z=[]
for atom in slab[id]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[id].translate([0.0,0.0,-shift])

mk_qe_input_relax(slab[id],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

view(slab[id])

os.chdir('../')



os.makedirs('./100',exist_ok=True)
os.chdir('./100')
slab = slabgen(cell,1,0,0,1,1,1,7.5,7.5)
#view(slab)
#sys.exit()
slab = slabgen(cell,1,0,0,nx100,ny100,nz100,7.5,7.5)
id = 0
#view(slab[0])
#view(slab[1])

c = FixAtoms(indices=[atom.index for atom in slab[0] if atom.position[2] < 9.5])
slab[0].set_constraint(c)

c = FixAtoms(indices=[atom.index for atom in slab[1] if atom.position[2] < 9.5])
slab[1].set_constraint(c)

#view(slab[0])
#view(slab[1])

z=[]
for atom in slab[0]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[0].translate([0.0,0.0,-shift])

z=[]
for atom in slab[1]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab[1].translate([0.0,0.0,-shift])

os.makedirs('surface0',exist_ok=True)
os.chdir('surface0')
mk_qe_input_relax(slab[0],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')
os.makedirs('surface1',exist_ok=True)
os.chdir('surface1')
mk_qe_input_relax(slab[1],'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

view(slab[0])
view(slab[1])

os.chdir('../')
# MoS2スラブモデル作成終了
"""

"""
# DOSプロット
#mpid='mp-1434'
mpid='mp-2815'
style='slab'
mill='001'
g.dos_rusult_dir = '/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/'+str(style)+'/dos/'+str(mpid)+'/'+str(mill)
os.makedirs('./MoS2_catalysis/'+str(style)+'/plot_dos/'+str(mpid)+'/'+str(mill),exist_ok=True)
os.chdir('./MoS2_catalysis/'+str(style)+'/plot_dos/'+str(mpid)+'/'+str(mill))
plot_qe_dos(g.dos_rusult_dir,-4.787,['Mo','S'],nspin=False)
# DOSプロット終了
"""

"""
#単体の構造最適化
ecutwfc0=77.0
ecutrho0=539.0
kpoint=4
scale=1.6

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

#mpid='S'
mpid='Mo'

cell = rd_cif(g.cifdir+'/'+str(mpid)+'.cif')

os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

os.makedirs('./bulk',exist_ok=True)
os.chdir('./bulk')

os.makedirs('./optimize',exist_ok=True)
os.chdir('./optimize')

os.makedirs('./'+str(mpid),exist_ok=True)
os.chdir('./'+str(mpid))

mk_qe_input_vcrelax(cell,'pbe','paw',level='high',nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
# 単体の構造最適化終了
"""

"""
# DOS計算
ecutwfc0=77.0
ecutrho0=539.0
kpoint=4
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

mpid = 1434
#mpid = 2815

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/bulk/optimize/mp'+str(mpid)+'/qe_vc-relax.pwo')

os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

os.makedirs('./bulk',exist_ok=True)
os.chdir('./bulk')

os.makedirs('./dos',exist_ok=True)
os.chdir('./dos')

os.makedirs('./mp'+str(mpid),exist_ok=True)
os.chdir('./mp'+str(mpid))

os.makedirs('./scf',exist_ok=True)
os.chdir('./scf')
mk_qe_input_scf(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
os.chdir('../')

os.makedirs('./dos',exist_ok=True)
os.chdir('./dos')
mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(12,12,12),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
# DOS計算終了
"""


"""
# バルクの構造最適化
ecutwfc0=77.0
ecutrho0=539.0
kpoints=[2,4,8]

scales=[1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0]

#mpid = 1434
mpid = 2815
cell = rd_cif(g.cifdir+'/MoS2_mp'+str(mpid)+'.cif')

os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

os.makedirs('./bulk',exist_ok=True)
os.chdir('./bulk')

os.makedirs('./mp'+str(mpid),exist_ok=True)
os.chdir('./mp'+str(mpid))

for scale in scales:
    os.makedirs('./scale_'+str(scale),exist_ok=True)
    os.chdir('./scale_'+str(scale))
    for k0 in kpoints:
        os.makedirs('./k_'+str(k0),exist_ok=True)
        os.chdir('./k_'+str(k0))
        ecutwfc=ecutwfc0*scale
        ecutrho=ecutrho0*scale
        os.makedirs('./opt',exist_ok=True)
        os.chdir('./opt')
        mk_qe_input_vcrelax(cell,'pbe','paw',level='high',nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
        os.chdir('../')
        #os.makedirs('./dos',exist_ok=True)
        #os.chdir('./dos')
        #mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(12,12,12),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
        #os.chdir('../')
        #os.makedirs('./band',exist_ok=True)
        #os.chdir('./band')
        #mk_qe_input_band(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
        #os.chdir('../')        
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
os.chdir('../')
os.chdir('../')
#バルクの構造最適化
"""

"""
# Pd原子の吸着モデル
# MoS2 + Pd ontop
os.makedirs('./MoS2_Pd1_ontop_test',exist_ok=True)
os.chdir('./MoS2_Pd1_ontop_test')
Pd = Atom('Pd',(3.192,3.686,4.909+2.5))
Pd1_ontop = slab001[id]+Pd
Pd1_ontop.write('Pd1_ontop.cif')
mk_qe_input_relax(Pd1_ontop,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pd bridge
os.makedirs('./MoS2_Pd1_bridge_test',exist_ok=True)
os.chdir('./MoS2_Pd1_bridge_test')
Pd = Atom('Pd',((3.192+6.384)/2.0,3.686,4.909+2.3))
Pd1_bridge = slab001[id]+Pd
Pd1_bridge.write('Pd1_bridge.cif')
mk_qe_input_relax(Pd1_bridge,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pd fcc
os.makedirs('./MoS2_Pd1_fcc_test',exist_ok=True)
os.chdir('./MoS2_Pd1_fcc_test')
Pd = Atom('Pd',(1.596,4.608,4.909+2.3))
Pd1_fcc = slab001[id]+Pd
Pd1_fcc.write('Pd1_fcc.cif')
mk_qe_input_relax(Pd1_fcc,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pd hcp
os.makedirs('./MoS2_Pd1_hcp_test',exist_ok=True)
os.chdir('./MoS2_Pd1_hcp_test')
Pd = Atom('Pd',((1.596+4.788+3.192)/3.0,(6.451+6.451+3.686)/3.0,4.909+2.3))
Pd1_hcp = slab001[id]+Pd
Pd1_hcp.write('Pd1_hcp.cif')
mk_qe_input_relax(Pd1_hcp,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pt ontop
os.makedirs('./MoS2_Pt1_ontop_test',exist_ok=True)
os.chdir('./MoS2_Pt1_ontop_test')
Pt = Atom('Pt',(3.192,3.686,4.909+2.5))
Pt1_ontop = slab001[id]+Pt
Pt1_ontop.write('Pt1_ontop.cif')
mk_qe_input_relax(Pt1_ontop,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pt bridge
os.makedirs('./MoS2_Pt1_bridge_test',exist_ok=True)
os.chdir('./MoS2_Pt1_bridge_test')
Pt = Atom('Pt',((3.192+6.384)/2.0,3.686,4.909+2.3))
Pt1_bridge = slab001[id]+Pt
Pt1_bridge.write('Pt1_bridge.cif')
mk_qe_input_relax(Pt1_bridge,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pt fcc
os.makedirs('./MoS2_Pt1_fcc_test',exist_ok=True)
os.chdir('./MoS2_Pt1_fcc_test')
Pt = Atom('Pt',(1.596,4.608,4.909+2.3))
Pt1_fcc = slab001[id]+Pt
Pt1_fcc.write('Pt1_fcc.cif')
mk_qe_input_relax(Pt1_fcc,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

# MoS2 + Pt hcp
os.makedirs('./MoS2_Pt1_hcp_test',exist_ok=True)
os.chdir('./MoS2_Pt1_hcp_test')
Pt = Atom('Pt',((1.596+4.788+3.192)/3.0,(6.451+6.451+3.686)/3.0,4.909+2.3))
Pt1_hcp = slab001[id]+Pt
Pt1_hcp.write('Pt1_hcp.cif')
mk_qe_input_relax(Pt1_hcp,'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')
# Pd原子の吸着モデル終了
"""

"""
# SIESTA
mpid = 1434
#mpid = 2815
cell = rd_cif(g.cifdir+'/MoS2_mp'+str(mpid)+'.cif')

os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

os.makedirs('./siesta',exist_ok=True)
os.chdir('./siesta')

slab001 = slabgen(cell,0,0,1,1,1,1,7.5,7.5)
id = 2

#view(slab001)
#view(slab001[id])

os.system('atomsk slab_3.cif -orthogonal-cell slab_3_ortho.cfg')
slab001 = read('slab_3_ortho.cfg')
slab001 = make_supercell(slab001,([5,0,0],[0,3,0],[0,0,1]),wrap=True)
lat = slab001.get_cell()
z=[]
for atom in slab001:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)
del slab001[[atom.index for atom in slab001 if atom.position[1] < 8.0 and atom.position[2] > 11.0]]
slab001.write('slab.xyz')
view(slab001)
smiles2xyz('Cl[Pd]Cl','PdCl2',False,userandom=True)
smiles2xyz('O','H2O',True)
smiles2xyz('CCO','CH3CH2OH',True)
smiles2xyz('[Cl-]','Cl',False)
smiles2xyz('[H+]','H',False)
mol_PdCl2 = read('PdCl2.xyz')
mol_H2O = read('H2O.xyz')
mol_CH3CH2OH = read('CH3CH2OH.xyz')
mol_Cl = read('Cl.xyz')
mol_H = read('H.xyz')
#view(mol_PdCl2)
#view(mol_H2O)
#view(mol_CH3CH2OH)
#view(mol_HCl)
buffer=2.0
slab_top=zmax+buffer
slab_bottom=zmin-buffer
molboxlist={
    'mol1':
        {
            'mol':'PdCl2',
            'num':1,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol2':
        {
            'mol':'H2O',
            'num':10,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol3':
        {
            'mol':'CH3CH2OH',
            'num':10,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol4':
        {
            'mol':'Cl',
            'num':2,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol5':
        {
            'mol':'H',
            'num':2,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[slab_top,lat[2][2]]
        },
    'mol6':
        {
            'mol':'PdCl2',
            'num':1,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol7':
        {
            'mol':'H2O',
            'num':10,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol8':
        {
            'mol':'CH3CH2OH',
            'num':10,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol9':
        {
            'mol':'Cl',
            'num':2,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        },
    'mol10':
        {
            'mol':'H',
            'num':2,
            'lx':[0.0,lat[0][0]],
            'ly':[0.0,lat[1][1]],
            'lz':[0.0,slab_bottom]
        }
}
mk_packmol_slab_random('slab',molboxlist,shift=0.5)
os.system('packmol < packmol.inp')

cell = read('system.xyz')
cell.set_cell(lat)

view(cell)
#print(len(cell))
#mk_qe_input_npt(cell,'pbe','paw',400,100,0,dt=0.5,level='high',estep=1000,nstep=200000,kpts=None,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
#mk_siesta_input_npt(cell,'VDW','DZP',100.0,[1,1,1],343.0,0.0,200000,pseudo_path=g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
mk_siesta_input_scf_withEfield(cell,'VDW','DZP',100.0,[1,1,1],pseudo_path=g.siesta_pot,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
"""