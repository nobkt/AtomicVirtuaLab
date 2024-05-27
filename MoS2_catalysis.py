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
from ase.optimize import BFGS
import numpy as np
import os
import sys

#os.system('export DISPLAY=localhost:11.0')

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir='./cifs'


from ocpmodels.common.relaxation.ase_utils import OCPCalculator
gnn='schnet'
def setocp():
    if gnn == 'dimnetpp':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/dimenetpp_all.pt"
        config_yml_path = "configs/s2ef/all/dimenet_plus_plus/dpp.yml"
        checkpoint_path = "content/ocp/dimenetpp_all.pt"
    elif gnn == 'gemnet-dT':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/gemnet_t_direct_h512_all.pt"
        config_yml_path = "configs/s2ef/all/gemnet/gemnet-dT.yml"
        checkpoint_path = "content/ocp/gemnet_t_direct_h512_all.pt"
    elif gnn == 'painn':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/painn_h512_s2ef_all.pt"
        config_yml_path = "configs/s2ef/all/painn/painn_h512.yml"
        checkpoint_path = "content/ocp/painn_h512_s2ef_all.pt"
    elif gnn == 'cgcnn':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/cgcnn_all.pt"
        config_yml_path = "configs/s2ef/all/cgcnn/cgcnn.yml"
        checkpoint_path = "content/ocp/cgcnn_all.pt"
    elif gnn == 'spinconv':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/spinconv_force_centric_all.pt"
        config_yml_path = "configs/s2ef/all/spinconv/spinconv_force.yml"
        checkpoint_path = "content/ocp/spinconv_force_centric_all.pt"
    elif gnn == 'schnet':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/schnet_all_large.pt"
        config_yml_path = "configs/s2ef/all/schnet/schnet.yml"
        checkpoint_path = "content/ocp/schnet_all_large.pt"
    elif gnn == 'gemnet_oc':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/gemnet_oc_base_s2ef_all_md.pt"
        config_yml_path = "configs/s2ef/all/gemnet/gemnet-oc.yml"
        checkpoint_path = "content/ocp/gemnet_oc_base_s2ef_all_md.pt"
    elif gnn == 'equiformer_v2':
        config_yml_file = "/home/modules/applications/common/ocp/configs/s2ef/all/*"
        checkpoint_file = "/home/modules/applications/common/ocp/content/ocp/eq2_153M_ec4_allmd.pt"
        config_yml_path = "configs/s2ef/all/equiformer_v2/equiformer_v2_N@20_L@6_M@3_153M.yml"
        checkpoint_path = "content/ocp/eq2_153M_ec4_allmd.pt"
    os.makedirs('./configs/s2ef/all',exist_ok=True)
    os.makedirs('./content/ocp',exist_ok=True)
    os.system('cp -r '+config_yml_file+' ./configs/s2ef/all')
    os.system('cp -r '+checkpoint_file+' ./content/ocp')
    ocpcalc = OCPCalculator(config_yml=config_yml_path, checkpoint_path=checkpoint_path)
    #atoms.calc=ocpcalc
    return ocpcalc


"""
# DOS calc
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

elm = 'Au'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = '001'

Pd_site = 'ontop'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

os.chdir(path)
Pd_site = 'fcc'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

os.chdir(path)
Pd_site = 'hcp'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)


os.chdir(path)
# mp-2815
mpid = 2815
mill = '001'

Pd_site = 'ontop'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

os.chdir(path)
Pd_site = 'fcc'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

os.chdir(path)
Pd_site = 'hcp'
slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm),exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/dos/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_'+str(elm))

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

mk_qe_input_scf(slab,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

# DOS calc 終了
"""


"""
# one Pd吸着後スラブ 吸着モデル ポテンシャル計算
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
Pd_site = 'fcc'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/Pd_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_potential',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_potential')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
rlist = [2.1,2.2,2.3,2.4,2.6,2.7,2.8,2.9]
adtype = 'Pd_ontop'
x0 = 3.154
y0 = 5.463
z0 = 9.018
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_1
rlist = [2.6,2.7,2.8,2.9,3.1,3.2,3.3,3.4,3.5]
adtype = 'S_ontop_1'
x0 = 4.828
y0 = 4.496
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_2
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'S_ontop_2'
x0 = -0.006
y0 = 7.287
z0 = 7.652
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_3
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'S_ontop_3'
x0 = -0.012
y0 = 1.808
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')


# bridge_1
rlist = [3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5]
adtype = 'bridge_1'
x0 = (3.154+1.480)/2.0
y0 = (7.396+4.496)/2.0
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge_2
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'bridge_2'
x0 = (3.154-0.006)/2.0
y0 = (7.396+7.287)/2.0
z0 = (7.707+7.652)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge_3
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'bridge_3'
x0 = (3.154+1.572)/2.0
y0 = (7.396+10.031)/2.0
z0 = (7.707+7.648)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# fcc_1
rlist = [1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9]
adtype = 'fcc_1'
x0 = (1.480-1.594-0.006)/3.0
y0 = (4.496+4.548+7.287)/3.0
z0 = (7.707+7.648+7.652)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# hcp_1
rlist = [2.6,2.7,2.8,2.9,3.1,3.2,3.3,3.4,3.5]
adtype = 'hcp_1'
x0 = (1.480+3.154-0.006)/3.0
y0 = (4.496+7.396+7.287)/3.0
z0 = (7.707+7.707+7.652)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# hcp_2
rlist = [1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9]
adtype = 'hcp_2'
x0 = (1.480-0.012-1.594)/3.0
y0 = (4.496+1.808+4.568)/3.0
z0 = (7.707+7.648+7.648)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

os.chdir(path)

# mp-2815
mpid = 2815
mill = '001'
Pd_site = 'hcp'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/Pd_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_potential',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_potential')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
rlist = [2.1,2.2,2.3,2.4,2.6,2.7,2.8,2.9]
adtype = 'Pd_ontop'
x0 = 1.576
y0 = 6.386
z0 = 12.102
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_1
rlist = [2.6,2.7,2.8,2.9,3.1,3.2,3.3,3.4,3.5]
adtype = 'S_ontop_1'
x0 = 3.249
y0 = 7.334
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_2
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'S_ontop_2'
x0 = 4.733
y0 = 4.545
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
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# S_ontop_3
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'S_ontop_3'
x0 = 3.157
y0 = 1.804
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')


# bridge_1
rlist = [3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5]
adtype = 'bridge_1'
x0 = (3.249+1.576)/2.0
y0 = (7.334+4.436)/2.0
z0 = (10.787+10.787)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge_2
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'bridge_2'
x0 = (1.576+4.733)/2.0
y0 = (4.436+4.545)/2.0
z0 = (10.787+10.736)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# bridge_3
rlist = [1.6,1.7,1.8,1.9,2.1,2.2,2.3,2.4]
adtype = 'bridge_3'
x0 = (1.576+3.157)/2.0
y0 = (4.436+1.804)/2.0
z0 = (10.787+10.732)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# fcc_1
rlist = [2.6,2.7,2.8,2.9,3.1,3.2,3.3,3.4,3.5]
adtype = 'fcc_1'
x0 = (1.576+3.249+4.733)/3.0
y0 = (4.436+7.334+4.545)/3.0
z0 = (10.787+10.787+10.736)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# fcc_2
rlist = [1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9]
adtype = 'fcc_2'
x0 = (1.576-0.005+3.157)/3.0
y0 = (4.436+1.804+1.084)/3.0
z0 = (10.787+10.732+10.732)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')

# hcp_1
rlist = [1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9]
adtype = 'hcp_1'
x0 = (3.249+4.733+6.319)/3.0
y0 = (7.334+4.545+7.281)/3.0
z0 = (10.787+10.736+10.732)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    os.makedirs('./r'+str(r),exist_ok=True)
    os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    os.chdir('../')
os.chdir('../')
# one Pd吸着後スラブ 吸着モデル ポテンシャル計算 終了
"""

"""
# one Pd吸着後スラブ 吸着モデル 構造最適化
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
Pd_site = 'fcc'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/Pd_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
rlist = [2.6]
adtype = 'Pd_ontop'
x0 = 3.154
y0 = 5.463
z0 = 9.018
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_1
rlist = [2.8]
adtype = 'S_ontop_1'
x0 = 4.828
y0 = 4.496
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_2
rlist = [2.2]
adtype = 'S_ontop_2'
x0 = -0.006
y0 = 7.287
z0 = 7.652
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_3
rlist = [2.2]
adtype = 'S_ontop_3'
x0 = -0.012
y0 = 1.808
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')


# bridge_1
rlist = [3.8]
adtype = 'bridge_1'
x0 = (3.154+1.480)/2.0
y0 = (7.396+4.496)/2.0
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# bridge_2
rlist = [1.8]
adtype = 'bridge_2'
x0 = (3.154-0.006)/2.0
y0 = (7.396+7.287)/2.0
z0 = (7.707+7.652)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# bridge_3
rlist = [1.7]
adtype = 'bridge_3'
x0 = (3.154+1.572)/2.0
y0 = (7.396+10.031)/2.0
z0 = (7.707+7.648)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_1
rlist = [1.5]
adtype = 'fcc_1'
x0 = (1.480-1.594-0.006)/3.0
y0 = (4.496+4.548+7.287)/3.0
z0 = (7.707+7.648+7.652)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_1
rlist = [3.1]
adtype = 'hcp_1'
x0 = (1.480+3.154-0.006)/3.0
y0 = (4.496+7.396+7.287)/3.0
z0 = (7.707+7.707+7.652)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_2
rlist = [1.6]
adtype = 'hcp_2'
x0 = (1.480-0.012-1.594)/3.0
y0 = (4.496+1.808+4.568)/3.0
z0 = (7.707+7.648+7.648)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

os.chdir(path)

# mp-2815
mpid = 2815
mill = '001'
Pd_site = 'hcp'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/'+str(mill)+'/Pd_optimize/'+str(Pd_site)+'/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
rlist = [2.6]
adtype = 'Pd_ontop'
x0 = 1.576
y0 = 6.386
z0 = 12.102
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_1
rlist = [2.8]
adtype = 'S_ontop_1'
x0 = 3.249
y0 = 7.334
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_2
rlist = [2.2]
adtype = 'S_ontop_2'
x0 = 4.733
y0 = 4.545
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# S_ontop_3
rlist = [2.2]
adtype = 'S_ontop_3'
x0 = 3.157
y0 = 1.804
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')


# bridge_1
rlist = [3.8]
adtype = 'bridge_1'
x0 = (3.249+1.576)/2.0
y0 = (7.334+4.436)/2.0
z0 = (10.787+10.787)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# bridge_2
rlist = [1.8]
adtype = 'bridge_2'
x0 = (1.576+4.733)/2.0
y0 = (4.436+4.545)/2.0
z0 = (10.787+10.736)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# bridge_3
rlist = [1.7]
adtype = 'bridge_3'
x0 = (1.576+3.157)/2.0
y0 = (4.436+1.804)/2.0
z0 = (10.787+10.732)/2.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_1
rlist = [3.1]
adtype = 'fcc_1'
x0 = (1.576+3.249+4.733)/3.0
y0 = (4.436+7.334+4.545)/3.0
z0 = (10.787+10.787+10.736)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_2
rlist = [1.6]
adtype = 'fcc_2'
x0 = (1.576-0.005+3.157)/3.0
y0 = (4.436+1.804+1.084)/3.0
z0 = (10.787+10.732+10.732)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_1
rlist = [1.5]
adtype = 'hcp_1'
x0 = (3.249+4.733+6.319)/3.0
y0 = (7.334+4.545+7.281)/3.0
z0 = (10.787+10.736+10.732)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')
# one Pd吸着後スラブ 吸着モデル 構造最適化 終了
"""

"""
###
# two Pd吸着後スラブ 吸着モデル 構造最適化
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
Pd_site = 'fcc'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd/Pd_optimize/fcc_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd2/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd2/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# fcc_1
rlist = [1.5]
adtype = 'fcc_1'
x0 = 6.316
y0 = 5.462
z0 = (7.648+7.708+7.652)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_2
rlist = [1.5]
adtype = 'fcc_2'
x0 = 4.735
y0 = 2.722
z0 = (7.640+7.708+7.644)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_3
rlist = [1.5]
adtype = 'fcc_3'
x0 = 4.739
y0 = 8.205
z0 = (7.695+7.648+7.641)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_4
rlist = [1.5]
adtype = 'fcc_4'
x0 = 1.575
y0 = 2.689
z0 = (7.698+7.638+7.640)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# fcc_5
rlist = [1.5]
adtype = 'fcc_5'
x0 = 1.576
y0 = 8.219
z0 = (7.694+7.695+7.645)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

os.chdir(path)

# mp-2815
mpid = 2815
mill = '001'
Pd_site = 'hcp'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_Pd/Pd_optimize/bridge_2/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd2/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd2/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# hcp_1
rlist = [1.5]
adtype = 'hcp_1'
x0 = -0.006
y0 = 9.120
z0 = (10.801+10.744+10.737)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_2
rlist = [1.5]
adtype = 'hcp_2'
x0 = -1.586
y0 = 6.384
z0 = (10.801+10.734+10.733)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_3
rlist = [1.5]
adtype = 'hcp_3'
x0 = 3.157
y0 = 9.126
z0 = (10.787+10.737+10.732)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_4
rlist = [1.5]
adtype = 'hcp_4'
x0 = -0.037
y0 = 3.634
z0 = (10.788+10.733+10.733)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# hcp_5
rlist = [1.5]
adtype = 'hcp_5'
x0 = 4.751
y0 = 6.397
z0 = (10.783+10.783+10.734)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')


# two Pd吸着後スラブ 吸着モデル 構造最適化 終了
###
"""

"""
# three Pd吸着後スラブ 吸着モデル 構造最適化
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
Pd_site = 'fcc'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd2/Pd_optimize/fcc_4/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd3/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd3/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# fcc_1
rlist = [1.5]
adtype = 'fcc_1'
x0 = 4.742
y0 = 2.723
z0 = (7.682+7.681+7.649)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# on_Pd3
rlist = [1.5]
adtype = 'on_Pd3'
x0 = 1.576
y0 = 4.548
z0 = (9.114+9.114+9.127)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = '001'
Pd_site = 'hcp'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_Pd2/Pd_optimize/hcp_4/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd3/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/on_'+str(Pd_site)+'_Pd3/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# hcp_1
rlist = [1.5]
adtype = 'hcp_1'
x0 = 4.741
y0 = 6.387
z0 = (10.772+10.772+10.740)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')

# on_Pd3
rlist = [1.5]
adtype = 'on_Pd3'
x0 = 1.577
y0 = 4.561
z0 = (12.214+12.218+12.215)/3.0
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',nosym=True,tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.chdir('../')


# tree Pd吸着後スラブ 吸着モデル 構造最適化 終了
###
"""

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

"""
# 吸着モデル ポテンシャル計算
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]

elm = 'H'

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
rlist = [1.1,1.2,1.3,1.4]
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
rlist = [0.6,0.7,0.8,0.9]
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
rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
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
rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
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
rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
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
rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
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
rlist = [1.0]
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
rlist = [1.0]
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
# 吸着モデル ポテンシャル計算終了
"""

"""
# 吸着モデル 構造最適化
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

#rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]

elm = 'H'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize_plusU',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize_plusU')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)

# ontop
#rlist = [2.2]
rlist = [1.3]
adtype = 'ontop'
x0 = 1.577
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# bridge
#rlist = [1.7]
rlist = [0.8]
adtype = 'bridge'
x0 = (1.577+4.731)/2.0
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc
#rlist = [1.5]
rlist = [1.0]
adtype = 'fcc'
x0 = 3.154
y0 = 5.463
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp
#rlist = [1.6]
rlist = [0.9]
adtype = 'hcp'
x0 = 1.577
y0 = 6.373
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    view(adsorp)
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize_plusU',exist_ok=True)
os.chdir('./MoS2_catalysis/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize_plusU')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# ontop
adtype = 'ontop'
#rlist = [2.2]
rlist = [1.4]
x0 = 1.576
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# bridge
#rlist = [1.7]
rlist = [0.8]
adtype = 'bridge'
x0 = (1.576+4.727)/2.0
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc
adtype = 'fcc'
#rlist = [1.6]
rlist = [0.9]
x0 = (1.576+4.727+3.151)/3.0
y0 = (4.548+4.548+7.278)/3.0
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp
#rlist = [1.5]
rlist = [1.0]
adtype = 'hcp'
x0 = 1.576
y0 = 6.368
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
for r in rlist:
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    #view(adsorp)
    #sys.exit()
    mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')
# 吸着モデル 構造最適化計算終了
"""


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
# DOSプロット 吸着モデル
#mpid='1434'
mpid='2815'
style='adsorp'
mill='001'
elm='H'
site='hcp'
g.dos_rusult_dir = '/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/dos/mp'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_adsorp/'+str(site)
os.makedirs('./MoS2_catalysis/'+str(style)+'/plot_dos/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_adsorp/'+str(site),exist_ok=True)
os.chdir('./MoS2_catalysis/'+str(style)+'/plot_dos/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_adsorp/'+str(site))
plot_qe_dos(g.dos_rusult_dir,-4.8992,['Mo','S','H'],nspin=False)
# DOSプロット 吸着モデル 終了
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
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

#mpid = 1434
mpid = 2815

adelm = 'H'
site = 'hcp_1'

#cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/bulk/optimize/mp'+str(mpid)+'/qe_vc-relax.pwo')
#cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/001/'+str(adelm)+'_optimize/'+str(site)+'/qe_relax.pwo')
#cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/001/on_fcc_Pd/'+str(adelm)+'_optimize/'+str(site)+'/qe_relax.pwo')
cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-'+str(mpid)+'/001/on_hcp_Pd/'+str(adelm)+'_optimize/'+str(site)+'/qe_relax.pwo')
view(cell)
os.makedirs('./MoS2_catalysis',exist_ok=True)
os.chdir('./MoS2_catalysis')

#os.makedirs('./bulk',exist_ok=True)
#os.chdir('./bulk')

os.makedirs('./adsorp',exist_ok=True)
os.chdir('./adsorp')

os.makedirs('./dos',exist_ok=True)
os.chdir('./dos')

os.makedirs('./mp'+str(mpid),exist_ok=True)
os.chdir('./mp'+str(mpid))

os.makedirs('./001',exist_ok=True)
os.chdir('./001')

#os.makedirs('./on_fcc_Pd',exist_ok=True)
#os.chdir('./on_fcc_Pd')
os.makedirs('./on_hcp_Pd',exist_ok=True)
os.chdir('./on_hcp_Pd')

os.makedirs('./'+str(adelm)+'_adsorp',exist_ok=True)
os.chdir('./'+str(adelm)+'_adsorp')

os.makedirs('./'+str(site),exist_ok=True)
os.chdir('./'+str(site))
mk_qe_input_scf(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,nosym=True,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
#mk_qe_input_scf(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
#mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,k0),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
os.chdir('../')

#os.makedirs('./dos',exist_ok=True)
#os.chdir('./dos')
#mk_qe_input_dos(cell,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(12,12,12),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4})
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

"""
# MoS2への吸着ポテンシャル OCP
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

potmin=[]
rlist=[]
for i in range(60):
    r = 0.5+float(i)*0.1
    if r <= 6.0:
        rlist.append(r)
#rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
elm = 'Pt'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)

# ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'ontop'
x0 = 1.577
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


# bridge
#rlist = [0.6,0.7,0.8,0.9]
adtype = 'bridge'
x0 = (1.577+4.731)/2.0
y0 = 4.552
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc
#rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
adtype = 'fcc'
x0 = 3.154
y0 = 5.463
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp
#rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
adtype = 'hcp'
x0 = 1.577
y0 = 6.373
z0 = 7.654
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = '001'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/slab/optimize/mp'+str(mpid)+'/001/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# ontop
#rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
adtype = 'ontop'
x0 = 1.576
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #sys.exit()
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# bridge
#rlist = [0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4]
adtype = 'bridge'
x0 = (1.576+4.727)/2.0
y0 = 4.548
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #sys.exit()
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc
#rlist = [1.0]
adtype = 'fcc'
x0 = (1.576+4.727+3.151)/3.0
y0 = (4.548+4.548+7.278)/3.0
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #sys.exit()
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp
#rlist = [1.0]
adtype = 'hcp'
x0 = 1.576
y0 = 6.368
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #sys.exit()
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir('../../..')
f = open(gnn+'_potmin.txt','a')
for d in potmin:
    f.write(d[0]+' '+d[1]+' '+str(d[2])+' '+str(d[3])+'\n')
f.close()
# MoS2への吸着ポテンシャル OCP 終了
"""

"""
# M-MoS2への吸着ポテンシャル
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

potmin=[]
rlist=[]
for i in range(60):
    r = 0.5+float(i)*0.1
    if r <= 6.0:
        rlist.append(r)
#rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
elm = 'Pt'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = 'Pd_on_fcc'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/Pd_optimize/fcc/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 3.154
y0 = 5.463
z0 = 9.018
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# near_hcp
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'near_hcp'
x0 = 1.577
y0 = 6.373
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = -0.001
y0 = 0.001
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = -1.578
y0 = 2.731
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = -3.154
y0 = 5.462
z0 = 7.645
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = -4.729
y0 = 8.195
z0 = 7.645
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 3.154
y0 = -0.002
z0 = 7.652
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_6'
x0 = 1.570
y0 = 2.718
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_7'
x0 = -0.015
y0 = 5.464
z0 = 7.652
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_8'
x0 = -1.579
y0 = 8.195
z0 = 7.645
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_9
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_9'
x0 = 6.309
y0 = 0.000
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_10
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_10'
x0 = 4.737
y0 = 2.718
z0 = 7.707
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_11
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_11'
x0 = 1.569
y0 = 8.206
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_12
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_12'
x0 = 9.461
y0 = 0.001
z0 = 7.643
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_13
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_13'
x0 = 7.886
y0 = 2.731
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_14
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_14'
x0 = 6.323
y0 = 5.464
z0 = 7.652
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_15
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_15'
x0 = 4.739
y0 = 8.206
z0 = 7.648
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = 'Pd_on_hcp'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/Pd_optimize/hcp/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 1.576
y0 = 6.366
z0 = 12.102
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# near_fcc
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'near_fcc'
x0 = (3.249+1.576+4.733)/3.0
y0 = (7.334+4.436+4.545)/3.0
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 1.576
y0 = 0.908
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hpc_2'
x0 = -0.009
y0 = 3.626
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = -1.591
y0 = 6.366
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = -3.153
y0 = 9.098
z0 = 10.729
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 4.726
y0 = 0.909
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_6'
x0 = 3.160
y0 = 3.626
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_7'
x0 = -0.006
y0 = 9.111
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_8'
x0 = 7.878
y0 = 0.912
z0 = 10.730
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_9
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_9'
x0 = 6.304
y0 = 3.638
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_10
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_10'
x0 = 4.743
y0 = 6.366
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_11
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_11'
x0 = 3.158
y0 = 9.111
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_12
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_12'
x0 = 11.030
y0 = 0.909
z0 = 10.730
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_13
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_13'
x0 = 9.452
y0 = 3.638
z0 = 10.730
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_14
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_14'
x0 = 7.878
y0 = 6.369
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_15
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_15'
x0 = 6.304
y0 = 9.098
z0 = 10.732
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')
"""

"""
# M2-MoS2への吸着ポテンシャル
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

potmin=[]
rlist=[]
for i in range(60):
    r = 0.5+float(i)*0.1
    if r <= 6.0:
        rlist.append(r)
#rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]
elm = 'Pt'

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = 'on_fcc_Pd2'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd/Pd_optimize/fcc_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 0.070
y0 = 5.492
z0 = 9.049
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pd_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_bridge'
x0 = (0.070+3.089)/2.0
y0 = (5.492+5.486)/2.0
z0 = 9.049
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = -1.586
y0 = 8.205
z0 = 7.694
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 1.576
y0 = 8.219
z0 = 7.695
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 6.316
y0 = 5.462
z0 = 7.708
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.735
y0 = 2.722
z0 = 7.708
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 1.575
y0 = 2.689
z0 = 7.698
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_PdPt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt/Pd_optimize/fcc_14/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 6.181
y0 = 5.499
z0 = 9.089
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 3.256
y0 = 5.486
z0 = 8.951
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# PdPt_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'PdPt_bridge'
x0 = (3.256+6.181)/2.0
y0 = (5.486+5.499)/2.0
z0 = 8.951
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 1.564
y0 = 8.208
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = -0.012
y0 = 5.463
z0 = 7.737
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 1.572
y0 = 2.721
z0 = 7.737
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.731
y0 = 2.679
z0 = 7.682
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 7.890
y0 = 2.726
z0 = 7.711
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_6'
x0 = -3.146
y0 = 5.461
z0 = 7.711
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_7'
x0 = -4.721
y0 = 8.207
z0 = 7.685
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_8'
x0 = 4.735
y0 = 8.223
z0 = 7.685
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_PdPt'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt/Pd_optimize/fcc_14/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 6.181
y0 = 5.499
z0 = 9.089
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 3.256
y0 = 5.486
z0 = 8.951
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# PdPt_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'PdPt_bridge'
x0 = (3.256+6.181)/2.0
y0 = (5.486+5.499)/2.0
z0 = 8.951
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 1.564
y0 = 8.208
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = -0.012
y0 = 5.463
z0 = 7.737
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 1.572
y0 = 2.721
z0 = 7.737
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.731
y0 = 2.679
z0 = 7.682
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 7.890
y0 = 2.726
z0 = 7.711
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_6'
x0 = -3.146
y0 = 5.461
z0 = 7.711
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_7'
x0 = -4.721
y0 = 8.207
z0 = 7.685
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_8'
x0 = 4.735
y0 = 8.223
z0 = 7.685
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pt2'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt/Pt_optimize/bridge_2/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 3.251
y0 = 5.565
z0 = 8.985
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pt_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_bridge'
x0 = (3.251+4.675)/2.0
y0 = (5.565+8.038)/2.0
z0 = 8.985
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 1.517
y0 = 8.210
z0 = 7.655
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = -0.009
y0 = 5.446
z0 = 7.744
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 1.571
y0 = 2.707
z0 = 7.744
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.737
y0 = 2.698
z0 = 7.705
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 6.334
y0 = 5.431
z0 = 7.705
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pt2'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt/Pt_optimize/fcc_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 1.728
y0 = 6.343
z0 = 12.067
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pt_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_bridge'
x0 = (1.728+4.583)/2.0
y0 = (6.343+6.345)/2.0
z0 = 12.067
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 3.152
y0 = 9.164
z0 = 10.736
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = -0.004
y0 = 9.103
z0 = 10.826
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = -1.584
y0 = 6.368
z0 = 10.826
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = -0.013
y0 = 3.625
z0 = 10.787
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 3.152
y0 = 3.608
z0 = 10.785
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_PdPt'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt/Pd_optimize/hcp_2/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 1.506
y0 = 6.291
z0 = 12.038
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 0.033
y0 = 3.769
z0 = 12.179
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# PtPd_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'PtPd_bridge'
x0 = (1.506+0.033)/2.0
y0 = (6.291+3.769)/2.0
z0 = 12.038
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 11.027
y0 = 0.902
z0 = 10.796
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 1.578
y0 = 0.903
z0 = 10.796
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 3.196
y0 = 3.613
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 4.738
y0 = 6.367
z0 = 10.821
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 3.157
y0 = 9.107
z0 = 10.821
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_6'
x0 = -0.006
y0 = 9.115
z0 = 10.797
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_7'
x0 = -1.603
y0 = 6.379
z0 = 10.797
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_8'
x0 = 9.438
y0 = 3.637
z0 = 10.769
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_PdPt'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt/Pd_optimize/hcp_2/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pt_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pt_ontop'
x0 = 1.506
y0 = 6.291
z0 = 12.038
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 0.033
y0 = 3.769
z0 = 12.179
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# PtPd_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'PtPd_bridge'
x0 = (1.506+0.033)/2.0
y0 = (6.291+3.769)/2.0
z0 = 12.038
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 11.027
y0 = 0.902
z0 = 10.796
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 1.578
y0 = 0.903
z0 = 10.796
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 3.196
y0 = 3.613
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 4.738
y0 = 6.367
z0 = 10.821
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 3.157
y0 = 9.107
z0 = 10.821
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_6
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_6'
x0 = -0.006
y0 = 9.115
z0 = 10.797
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_7
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_7'
x0 = -1.603
y0 = 6.379
z0 = 10.797
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_8
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_8'
x0 = 9.438
y0 = 3.637
z0 = 10.769
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pd2'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_Pd/Pd_optimize/bridge_2/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# Pd_ontop
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_ontop'
x0 = 1.637
y0 = 6.339
z0 = 12.141
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# Pd_bridge
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'Pd_bridge'
x0 = (1.637+3.138)/2.0
y0 = (6.339+3.731)/2.0
z0 = 12.141
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = -0.037
y0 = 3.634
z0 = 10.788
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = -1.586
y0 = 6.384
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = -0.006
y0 = 9.120
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 3.157
y0 = 9.126
z0 = 10.783
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 4.751
y0 = 6.397
z0 = 10.783
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')
"""




# M3-MoS2への吸着ポテンシャル
ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

potmin=[]
rlist=[]
for i in range(60):
    r = 0.5+float(i)*0.1
    if r <= 6.0:
        rlist.append(r)
#rlist = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0]

path = os.getcwd()
os.chdir(path)

# mp-1434
mpid = 1434
mill = 'on_fcc_Pt3'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt2/Pt_optimize/fcc_11/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt3'
x0 = 3.126
y0 = 7.285
z0 = 9.031
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 6.318
y0 = 5.442
z0 = 7.717
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 4.737
y0 = 2.702
z0 = 7.717
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pd3'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd2/Pd_optimize/fcc_4/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd3'
x0 = 1.576
y0 = 4.548
z0 = 9.114
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 4.742
y0 = 2.723
z0 = 9.114
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 6.320
y0 = 5.459
z0 = 9.114
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pd2Pt'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd2/Pt_optimize/fcc_5/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd2Pt
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd2Pt'
x0 = 1.576
y0 = 4.609
z0 = 9.115
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 3.160
y0 = -0.012
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 4.745
y0 = 2.724
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 6.321
y0 = 5.460
z0 = 7.671
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.734
y0 = 8.205
z0 = 7.680
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 1.576
y0 = 8.205
z0 = 7.680
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pd2Pt'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pd2/Pt_optimize/fcc_5/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd2Pt
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd2Pt'
x0 = 1.576
y0 = 4.609
z0 = 9.115
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 3.160
y0 = -0.012
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 4.745
y0 = 2.724
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 6.321
y0 = 5.460
z0 = 7.671
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = 4.734
y0 = 8.205
z0 = 7.680
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = 1.576
y0 = 8.205
z0 = 7.680
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pt2Pd'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt2/Pd_optimize/fcc_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt2Pd
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt2Pd'
x0 = 3.091
y0 = 7.305
z0 = 9.027
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 6.321
y0 = 5.440
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 4.737
y0 = 2.702
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 1.565
y0 = 2.707
z0 = 7.710
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = -0.015
y0 = 5.445
z0 = 7.722
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = -1.593
y0 = 8.178
z0 = 7.660
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-1434
mpid = 1434
mill = 'on_fcc_Pt2Pd'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-1434/001/on_fcc_Pt2/Pd_optimize/fcc_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = -3.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt2Pd
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt2Pd'
x0 = 3.091
y0 = 7.305
z0 = 9.027
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_1'
x0 = 6.321
y0 = 5.440
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_2'
x0 = 4.737
y0 = 2.702
z0 = 7.712
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_3'
x0 = 1.565
y0 = 2.707
z0 = 7.710
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_4'
x0 = -0.015
y0 = 5.445
z0 = 7.722
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# fcc_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'fcc_5'
x0 = -1.593
y0 = 8.178
z0 = 7.660
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pt3'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt2/Pt_optimize/hcp_10/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt3'
x0 = 3.152
y0 = 7.310
z0 = 12.117
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 3.152
y0 = 3.627
z0 = 10.800
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = -0.008
y0 = 3.626
z0 = 10.800
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pd3'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_Pd2/Pd_optimize/hcp_4/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd3'
x0 = 1.577
y0 = 4.561
z0 = 12.215
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 4.741
y0 = 6.387
z0 = 10.722
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 6.313
y0 = 3.651
z0 = 10.722
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pt2Pd'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt2/Pd_optimize/hcp_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt2Pd
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt2Pd'
x0 = 3.152
y0 = 7.349
z0 = 12.303
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 3.152
y0 = 3.624
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 6.312
y0 = 3.627
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 7.892
y0 = 6.373
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 6.314
y0 = 9.109
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 7.870
y0 = 0.922
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pt2Pd'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hpc_Pt2/Pd_optimize/hcp_1/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pt2Pd
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pt2Pd'
x0 = 3.152
y0 = 7.349
z0 = 12.303
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 3.152
y0 = 3.624
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 6.312
y0 = 3.627
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 7.892
y0 = 6.373
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 6.314
y0 = 9.109
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 7.870
y0 = 0.922
z0 = 10.793
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pd2Pt'
elm = 'Pd'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_PdPt/Pd_optimize/hcp_3/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd2Pt
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd2Pt'
x0 = 1.583
y0 = 4.489
z0 = 12.096
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 1.577
y0 = 0.897
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 4.731
y0 = 0.898
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 6.315
y0 = 3.640
z0 = 10.757
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 4.744
y0 = 6.373
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 3.160
y0 = 9.107
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')


os.chdir(path)
# mp-2815
mpid = 2815
mill = 'on_hcp_Pd2Pt'
elm = 'Pt'

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/MoS2_catalysis/absorp/optimize/mp-2815/001/on_hcp_PdPt/Pd_optimize/hcp_3/qe_relax.pwo')

os.makedirs('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize',exist_ok=True)
os.chdir('./MoS2_catalysis_ocp/adsorp/mp-'+str(mpid)+'/'+str(mill)+'/'+str(elm)+'_optimize')

#view(slab)
lowpos = 0.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

#view(slab)
#sys.exit()

# on_Pd2Pt
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'on_Pd2Pt'
x0 = 1.583
y0 = 4.489
z0 = 12.096
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_1
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_1'
x0 = 1.577
y0 = 0.897
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_2
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_2'
x0 = 4.731
y0 = 0.898
z0 = 10.766
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_3
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_3'
x0 = 6.315
y0 = 3.640
z0 = 10.757
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_4
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_4'
x0 = 4.744
y0 = 6.373
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')

# hcp_5
#rlist = [1.1,1.2,1.3,1.4]
adtype = 'hcp_5'
x0 = 3.160
y0 = 9.107
z0 = 10.801
os.makedirs('./'+str(adtype),exist_ok=True)
os.chdir('./'+str(adtype))
ocpcalc = setocp()
pelist=[]
for r in rlist:
    #os.makedirs('./r'+str(r),exist_ok=True)
    #os.chdir('./r'+str(r))
    atm = Atom(elm,(x0,y0,z0+float(r)))
    adsorp = slab+atm
    adsorp.calc = ocpcalc
    pe = adsorp.get_potential_energy()
    f = open(gnn+'_potential.txt','a')
    f.write(str(r)+' '+str(pe)+'\n')
    f.close()
    pelist.append([r,pe])
    #view(adsorp)
    #mk_qe_input_scf(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #mk_qe_input_dos(adsorp,'pbe','paw',level='high',ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
    #os.chdir('../')
os.system('rm -rf ./configs')
os.system('rm -rf ./content')
pelist.sort(reverse=False, key=lambda x:x[1])
potmin.append(['mp'+str(mpid),adtype,pelist[0][0],pelist[0][1]])
atm = Atom(elm,(x0,y0,z0+float(pelist[0][0])))
adsorp = slab+atm
view(adsorp)
mk_qe_input_relax(adsorp,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)
os.chdir('../')






