from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.io import rd_cif, rd_lammpsdata
from AtomicVirtuaLab.deepmd import qe2dp, get_deepmd_list, wt_deepmd_json
from AtomicVirtuaLab.lammps import mk_npt_melt_input_deepmd, mk_nvt_input_uff,mk_npt_input_deepmd
from ase.build import make_supercell
from ase.visualize import view
from ase.units import mol
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.tools import add_displacement, deform_cell, scale_cell
import os
import shutil
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"


"""
# 学習データ作成 GPAW
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')

os.makedirs('a-Al2O3',exist_ok=True)
os.chdir('a-Al2O3')

ldens = 2.93

nx=2
ny=2
nz=1

cell = rd_cif(g.cifs+'/'+'Al2O3.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)

lat = cell.get_cell_lengths_and_angles()
volume = cell.get_volume()
masses = cell.get_masses()
mass = sum(masses)
density = (mass/mol)/(volume*1.0e-24)
scale = density/ldens
scale = (scale)**(1/3)

lat[0] = lat[0]*scale
lat[1] = lat[1]*scale
lat[2] = lat[2]*scale

P=[0.0,50000.0,-50000.0]
T=[300,1000,2000,3000,4000,5000]


for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        #view(cell_solid)
        #sys.exit()
        cell_liquid = cell.copy()
        cell_liquid.set_cell(lat,scale_atoms=True)
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #view(cell_liquid)
        #sys.exit()
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        #mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,100.0,12345)
        os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        view(cell_liquid)
        os.chdir('../')
        os.system('rm -rf tmp')
        #mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
# 学習データ作成 GPAW 終了
"""

"""
# 学習データ作成
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')

os.makedirs('a-Al2O3',exist_ok=True)
os.chdir('a-Al2O3')

ldens = 2.93

nx=2
ny=2
nz=1

cell = rd_cif(g.cifs+'/'+'Al2O3.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)

lat = cell.get_cell_lengths_and_angles()
volume = cell.get_volume()
masses = cell.get_masses()
mass = sum(masses)
density = (mass/mol)/(volume*1.0e-24)
scale = density/ldens
scale = (scale)**(1/3)

lat[0] = lat[0]*scale
lat[1] = lat[1]*scale
lat[2] = lat[2]*scale

P=[0]
T=[100,500,1000,2000,3000,4000,5000]

for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        #view(cell_solid)
        #sys.exit()
        cell_liquid = cell.copy()
        cell_liquid.set_cell(lat,scale_atoms=True)
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #view(cell_liquid)
        #sys.exit()
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,100.0,12345)
        os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        view(cell_liquid)
        os.chdir('../')
        os.system('rm -rf tmp')
        mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=2000,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
# 学習データ作成 終了
"""

"""
# 学習データ変換
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')
os.makedirs('train',exist_ok=True)
os.chdir('./train')
datadir = '/home/A23321P/work/myDeepMD/Al2O3/datas'
qe2dp(datadir,'pwo')
# 学習データ変換 終了
"""

"""
# ポテンシャル学習
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')
os.makedirs('train',exist_ok=True)
os.chdir('./train')
dpdir = './deepmd'
dp_list=get_deepmd_list(dpdir)
wt_deepmd_json(dpdir,dp_list,8.0,1000000,prec='high')
# ポテンシャル学習 終了
"""

"""
# 融点計算
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')
os.makedirs('./md',exist_ok=True)
os.chdir('./md')
os.makedirs('./melt',exist_ok=True)
os.chdir('./melt')

cell = rd_cif(g.cifs+'/'+'Al2O3.cif')
cell = make_supercell(cell,([5,0,0],[0,5,0],[0,0,4]),wrap=True)
lat = cell.get_cell()
z0 = lat[2][2]/2.0
for T0 in [2310,2320,2330,2340,2350,2360,2370,2380,2390]:
    os.makedirs('./T'+str(T0),exist_ok=True)
    os.chdir('./T'+str(T0))
    shutil.copy(g.forcedir+'/Al2O3_graph.pb','./graph.pb')
    mk_npt_melt_input_deepmd(cell,0.0005,10000,10000,200000,200000,10000000,1000,5000,T0,0.0,z0,12345)
    os.chdir('../')
# 融点計算終了
"""

# 速度検証
os.makedirs('./Al2O3_growth',exist_ok=True)
os.chdir('./Al2O3_growth')
os.makedirs('./benchmark',exist_ok=True)
os.chdir('./benchmark')
cell = rd_cif(g.cifs+'/'+'Al2O3.cif')
cell = make_supercell(cell,([8,0,0],[0,8,0],[0,0,4]),wrap=True)
print(len(cell))
view(cell)
shutil.copy(g.forcedir+'/Al2O3_graph.pb','./graph.pb')
mk_npt_input_deepmd(cell,0.0005,1,1,10,300,1.0,12345,mol=False)
#mk_npt_melt_input_deepmd(cell,0.0005,10000,10000,200000,200000,10000000,1000,5000,T0,0.0,z0,12345)
# 速度検証終了