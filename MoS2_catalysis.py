import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.io import rd_cif, smiles2xyz
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.packmol import mk_packmol_slab_random
from AtomicVirtuaLab.espresso import mk_qe_input_npt, mk_qe_input_relax, mk_qe_input_vcrelax, mk_qe_input_dos, mk_qe_input_band,mk_qe_input_scf,plot_qe_dos
from AtomicVirtuaLab.siesta import mk_siesta_input_npt
from ase.io import read, write
from ase.build import make_supercell, add_adsorbate
from ase.visualize import view
from ase import Atom
from ase.constraints import FixAtoms
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'

g.cifdir='./cifs'

# plot dos
#mpid='mp1434'
mpid='mp2815'
g.dos_rusult_dir = '/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/bulk/dos/'+str(mpid)
os.makedirs('./MoS2_catalysis/bulk/plot_dos/'+str(mpid),exist_ok=True)
os.chdir('./MoS2_catalysis/bulk/plot_dos/'+str(mpid))
plot_qe_dos(g.dos_rusult_dir,11.8780,['Mo','S'],nspin=False)



"""
#単体
# DOS
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
"""

"""
# DOS
ecutwfc0=77.0
ecutrho0=539.0
kpoint=4
scale=1.6

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
"""


"""
# bulk
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
"""

"""
# <001> surface
os.makedirs('./001',exist_ok=True)
os.chdir('./001')
slab001 = slabgen(cell,0,0,1,4,4,1,7.5,7.5)
id = 0

c = FixAtoms(indices=[atom.index for atom in slab001[id] if atom.position[2] < 12])
slab001[id].set_constraint(c)
view(slab001[id])

z=[]
for atom in slab001[id]:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab001[id].translate([0.0,0.0,-shift])

# MoS2 001
os.makedirs('./MoS2_001_test',exist_ok=True)
os.chdir('./MoS2_001_test')
view(slab001[id])
slab001[id].write('slab001.cif')
mk_qe_input_relax(slab001[id],'pbe','paw',level='high',ecutwfc=77.0,ecutrho=539.0,kpts=None,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nosym=True)
os.chdir('../')

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
slab001 = read('slab_1_ortho.cfg')
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
mk_siesta_input_npt(cell,'VDW','DZP',100.0,[1,1,1],343.0,0.0,200000,pseudo_path=g.siesta_pot,SolutionMethod='diagon',MixingWeight=0.1,MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
"""