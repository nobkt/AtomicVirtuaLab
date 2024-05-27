import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.io import rd_cif
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.espresso import mk_qe_input_vcrelax, mk_qe_input_relax
from AtomicVirtuaLab.lammps import mk_npt_input_deepmd
from ase.io import read
from ase.visualize import view
from ase.constraints import FixAtoms
from ase.build import make_supercell
from ase.cluster.cubic import FaceCenteredCubic
from ase import Atom
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
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



# M-Mo2C吸着 スキャン
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

mill = '111-MoC_1'
elm = 'Pt'
ngrid = 100
z0 = 1.283

path = os.getcwd()
os.chdir(path)

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/Mo2C_catalysis/slab/'+str(mill)+'/qe_relax.pwo')

#view(slab)
#sys.exit()

os.makedirs('./Mo2C_catalysis/adsorp/'+str(mill)+'/'+str(elm)+'_scan',exist_ok=True)
os.chdir('./Mo2C_catalysis/adsorp/'+str(mill)+'/'+str(elm)+'_scan')

#view(slab)
lowpos = -2.0

c = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < lowpos])
slab.set_constraint(c)

lat = slab.get_cell()

#view(slab)
#sys.exit()

nstep = 0
for xi in range(ngrid):
    for yi in range(ngrid):
        x0 = lat[0][0]*float(xi)/float(ngrid)+lat[0][1]*float(yi)/float(ngrid)
        y0 = lat[1][0]*float(xi)/float(ngrid)+lat[1][1]*float(yi)/float(ngrid)
        ocpcalc = setocp()
        pelist=[]
        for r in rlist:
            #os.makedirs('./r'+str(r),exist_ok=True)
            #os.chdir('./r'+str(r))
            atm = Atom(elm,(x0,y0,z0+float(r)))
            adsorp = slab+atm
            adsorp.calc = ocpcalc
            pe = adsorp.get_potential_energy()
            pelist.append([r,pe])
            print('step:',nstep,x0,y0,r,pe)
            nstep = nstep + 1
        os.system('rm -rf ./configs')
        os.system('rm -rf ./content')
        pelist.sort(reverse=False, key=lambda x:x[1])
        potmin.append([x0,y0,pelist[0][0],pelist[0][1]])

f = open(gnn+'_'+mill+'_'+elm+'_potmin.txt','a')
for d in potmin:
    f.write(str(d[0])+' '+str(d[1])+' '+str(d[2])+' '+str(d[3])+'\n')
f.close()

# M-Mo2C吸着終了 スキャン


"""
# OC20テスト
os.makedirs('./Mo2C_catalysis',exist_ok=True)
os.chdir('./Mo2C_catalysis')
os.makedirs('OC20test',exist_ok=True)
os.chdir('./')

slab = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/Mo2C_catalysis/slab/002/qe_relax.pwo')
slab = make_supercell(slab,([4,0,0],[0,3,0],[0,0,1]),wrap=False)
lat = slab.get_cell()
lat[2][2] = 50.0
slab.set_cell(lat)
z=[]
for atom in slab:
    z.append(atom.position[2])
zmax = max(z)
com1 = slab.get_center_of_mass()
view(slab)

surfaces = [(1, 0, 0), (1, 1, 1), (1, -1, 1)]
layers = [6, 5, -1]
trunc = FaceCenteredCubic('Pt', surfaces, layers)
trunc.rotate(45, 'x', rotate_cell=True)
trunc.rotate(54.736, 'z', rotate_cell=True)
trunc.rotate(90, 'y', rotate_cell=True)
z=[]
for atom in trunc:
    z.append(atom.position[2])
zmin = min(z)
com2 = trunc.get_center_of_mass()
view(trunc)

shift=[com1[0]-com2[0],com1[1]-com2[1],zmax+2.5-zmin]
trunc.translate(shift)

adsorp = slab+trunc
view(adsorp)

mk_npt_input_deepmd(adsorp,0.0005,10,10,200000,300,0.0,12345,mol=False)

# OC20テスト 終了
"""

"""
# スラブの構造最適化

ecutwfc0=77.0
ecutrho0=539.0
kpoint=2
scale=1.0

ecutwfc=77.0*scale
ecutrho=539.0*scale
k0 = kpoint

mpid = 1552

#110-Mo
nx110_Mo = 3
ny110_Mo = 2
nz110_Mo = 2
zlow_110_Mo = 11.0

#110-Mo/C_2
nx110_MoC_2 = 3
ny110_MoC_2 = 2
nz110_MoC_2 = 2
zlow_110_MoC_2 = 11.0

#100-Mo
nx100_Mo = 2
ny100_Mo = 2
nz100_Mo = 2
zlow_100_Mo = 12.0

#111-Mo/C_1
nx111_MoC_1 = 2
ny111_MoC_1 = 2
nz111_MoC_1 = 2
zlow_111_MoC_1 = 10.5

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/pwos/Mo2C_catalysis/bulk/mp'+str(mpid)+'/opt/qe_vc-relax.pwo')

os.makedirs('./Mo2C_catalysis/slab/mp-'+str(mpid)+'/opt',exist_ok=True)
os.chdir('./Mo2C_catalysis/slab/mp-'+str(mpid)+'/opt')

#slab = slabgen(cell,1,1,1,2,2,1,7.5,7.5)
#view(slab)
#sys.exit()

#110-Mo
os.makedirs('./110-Mo',exist_ok=True)
os.chdir('./110-Mo')

slab = slabgen(cell,1,1,0,nx110_Mo,ny110_Mo,nz110_Mo,7.5,7.5)
slab_ = slab[2].copy()

#print(len(slab_))
#view(slab_)
#sys.exit()

c = FixAtoms(indices=[atom.index for atom in slab_ if atom.position[2] < zlow_110_Mo])
slab_.set_constraint(c)
#view(slab_)

z=[]
for atom in slab_:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab_.translate([0.0,0.0,-shift])

mk_qe_input_relax(slab_,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

#view(slab_)

os.chdir('../')

#110-Mo/C_2
os.makedirs('./110-MoC_2',exist_ok=True)
os.chdir('./110-MoC_2')

slab = slabgen(cell,1,1,0,nx110_MoC_2,ny110_MoC_2,nz110_MoC_2,7.5,7.5)
slab_ = slab[3].copy()

#print(len(slab_))
#view(slab_)
#sys.exit()

c = FixAtoms(indices=[atom.index for atom in slab_ if atom.position[2] < zlow_110_MoC_2])
slab_.set_constraint(c)
#view(slab_)

z=[]
for atom in slab_:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab_.translate([0.0,0.0,-shift])

mk_qe_input_relax(slab_,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

#view(slab_)

os.chdir('../')

#100-Mo
os.makedirs('./100-Mo',exist_ok=True)
os.chdir('./100-Mo')

slab = slabgen(cell,1,0,0,nx100_Mo,ny100_Mo,nz100_Mo,7.5,7.5)
slab_ = slab[1].copy()

#print(len(slab_))
#view(slab_)
#sys.exit()

c = FixAtoms(indices=[atom.index for atom in slab_ if atom.position[2] < zlow_100_Mo])
slab_.set_constraint(c)
#view(slab_)

z=[]
for atom in slab_:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab_.translate([0.0,0.0,-shift])

mk_qe_input_relax(slab_,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

#view(slab_)

os.chdir('../')

#111-Mo/C_1
os.makedirs('./111-MoC_1',exist_ok=True)
os.chdir('./111-MoC_1')

slab = slabgen(cell,1,1,1,nx111_MoC_1,ny111_MoC_1,nz111_MoC_1,7.5,7.5)
slab_ = slab[3].copy()

#print(len(slab_))
#view(slab_)
#sys.exit()

c = FixAtoms(indices=[atom.index for atom in slab_ if atom.position[2] < zlow_111_MoC_1])
slab_.set_constraint(c)
#view(slab_)

z=[]
for atom in slab_:
    z.append(atom.position[2])
zmax=max(z)
zmin=min(z)

shift = (zmax + zmin)/2.0

slab_.translate([0.0,0.0,-shift])

mk_qe_input_relax(slab_,'pbe','paw',level='high',estep=1000,nstep=1000,nosym=False,ecutwfc=ecutwfc,ecutrho=ecutrho,mixing_beta=0.2,kpts=(k0,k0,1),ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4,'assume_isolated':'esm','esm_bc':'bc1'},nspin=False)

#view(slab_)

os.chdir('../')

# スラブの構造最適化 終了
"""


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