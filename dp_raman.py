from ase.build import molecule
from ase.io import read, write
from ase.visualize import view
from ase.build import sort, make_supercell
from AtomicVirtuaLab.espresso import mk_qe_input_relax, mk_qe_input_npt
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize,mk_siesta_input_nvt,mk_siesta_input_scf_withEfield,mk_siesta_input_scf_withEfield_wannier,mk_siesta_input_cellopt,get_valence,mk_siesta_input_npt
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_nvt_input_uff_rigid_scale, mk_npt_input_deepmd
from AtomicVirtuaLab.deepmd import get_deepmd_list, wt_deepmd_json, siesta2dp
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol
from AtomicVirtuaLab.io import rd_lammpsdata_init
import shutil
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"

"""
# DeepMD NPT
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

os.makedirs('./deepmd',exist_ok=True)
os.chdir('./deepmd')

#cell = rd_lammpsdata_init({1:35,2:6,3:7,4:30},g.cifdir+'/Zn-Pc-allBr_64.data',True)
cell = read(g.cifdir+'/Zn-Pc-allBr_64.cif')
view(cell)
cell = sortmol(cell)

shutil.copy(g.forcedir+'/ZnPcBr_SZ_graph.pb','./graph.pb')
mk_npt_input_deepmd(cell,0.0005,2000,2000,2000000,300,0.0,12345,mol=True)


# DeepMD NPT 終了
"""

"""
# ポテンシャル学習
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('train',exist_ok=True)
os.chdir('./train')
dpdata=True
if dpdata:
    path = os.getcwd()
    datadir_ = '/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-allBr/datas/KBM'
    os.chdir(datadir_)
    siesta2dp()
    os.system('cp -r '+str(datadir_)+'/deepmd '+str(path))
    sys.exit()
dpdir = './deepmd'
dp_list=get_deepmd_list(dpdir)
wt_deepmd_json(dpdir,dp_list,8.0,1000000,prec='high')
# ポテンシャル学習 終了
"""

"""
# dipleおよびpolarizability解析
#test
C = 1.602176634e-19
A2m = 1.0e-10
Cm2Deb = 3.33564e-30
dir_ = '/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-allBr/wannier/DZP/cutoff_100.0/ez-'

atomic_mu=[]
grobal_mu=[]

wxyz = read(dir_+'/siesta_centres.xyz')
strout = read(dir_+'/siesta.STRUCT_OUT',format='struct_out')
lat = strout.get_cell()

wxyz.set_cell(lat)

wxyz = sortmol(wxyz)
nmol = max(wxyz.arrays['mol-id'])
wxyz.set_tags(wxyz.arrays['mol-id'])

wxyz.write('test.cif')

valence = get_valence(dir_)

mu1=[]
grobal_mux = 0.0
grobal_muy = 0.0
grobal_muz = 0.0
imol = 0
for imol0 in range(nmol):
    imol = imol0 + 1
    ux = 0.0
    uy = 0.0
    uz = 0.0
    for atom in wxyz:
        if atom.tag == imol:
            if atom.symbol != 'X':
                ux = ux + valence[atom.symbol]*atom.position[0]*C*A2m/Cm2Deb
                uy = uy + valence[atom.symbol]*atom.position[1]*C*A2m/Cm2Deb
                uz = uz + valence[atom.symbol]*atom.position[2]*C*A2m/Cm2Deb
            elif atom.symbol == 'X':
                ux = ux - 2.0*atom.position[0]*C*A2m/Cm2Deb
                uy = uy - 2.0*atom.position[1]*C*A2m/Cm2Deb
                uz = uz - 2.0*atom.position[2]*C*A2m/Cm2Deb
    mu1.append([ux,uy,uz])
    grobal_mux = grobal_mux + ux
    grobal_muy = grobal_muy + uy
    grobal_muz = grobal_muz + uz

grobal_mu.append([grobal_mux,grobal_muy,grobal_muz])

atomic_mu.append(mu1)

print(atomic_mu[0])
print(grobal_mu[0])
print(sum(grobal_mu[0]))

# dipleおよびpolarizability解析終了
"""

"""
# Zn-Pc アモルファスWannier
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('Zn-Pc-amo',exist_ok=True)
os.chdir('Zn-Pc-amo')
os.makedirs('wannier',exist_ok=True)
os.chdir('wannier')

#cell = read(g.cifdir+'/Zn-Pc_10.data',format='lammps-data',sort_by_id=True,Z_of_type={1:6,2:1,3:7,4:30})
#cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-amo/8/cell-opt/siesta.STRUCT_OUT',format='struct_out')
cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-allBr-amo/cell-opt/siesta.STRUCT_OUT',format='struct_out')

view(cell)

#mk_siesta_input_optimize(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

#mk_siesta_input_cellopt(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

for bandscale in range(31):
    os.makedirs('./'+str(bandscale),exist_ok=True)
    os.chdir('./'+str(bandscale))
    mk_siesta_input_scf_withEfield_wannier(cell,'PBE','SZ',50.0,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')


# Zn-Pc アモルファスWannier 終了
"""

"""
# Zn-Pc アモルファスモデル読み込み
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('Zn-Pc-allBr-amo',exist_ok=True)
os.chdir('Zn-Pc-allBr-amo')

cell = read(g.cifdir+'/Zn-Pc-allBr_8.data',format='lammps-data',sort_by_id=True,Z_of_type={1:35,2:6,3:7,4:30})

view(cell)

#mk_siesta_input_optimize(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

#mk_siesta_input_cellopt(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

# Zn-Pc アモルファスモデル読み込み 終了
"""

"""
# Zn-Pc 結晶モデル作成
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./cif_test',exist_ok=True)
os.chdir('./cif_test')

cell = read(g.cifdir+'/ZincPhthalocyanine.cif')
cell_solid = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/cell_opt/Zn-Pc_unitcell/siesta.STRUCT_OUT',format='struct_out')
cell_test = read('POSCAR')
view(cell_test)
#view(cell_solid)
cell_solid = make_supercell(cell,([2,0,0],[0,2,0],[0,0,2]),wrap=True)
#view(cell_solid)

#print(len(cell))
#cell = read(g.cifdir+'/Zn-Pc_2mer.xyz')
#cell.set_cell([50,50,50])

cell = sortmol(cell)
cell_solid = sortmol(cell_solid)
cell_test = sortmol(cell_test)

#mk_siesta_input_cellopt(cell,'PBE','SZ',50.0,[1,2,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
mk_siesta_input_scf_withEfield_wannier(cell_solid,'PBE','SZ',50.0,None,g.siesta_pot,bandscale=30,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

#cell.write('test.lmp',format='lammps-data',atom_style='full')
#view(cell)

# Zn-Pc 結晶モデル作成 終了
"""


"""
# 双極子モーメントと分極率テンソルの試計算
au2v=51.4220632
epsilon=0.001*au2v
xc='KBM'
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./dipole_and_polarizability',exist_ok=True)
os.chdir('./dipole_and_polarizability')

path = os.getcwd()
os.chdir(path)

# Zn-Pc-allBr
os.chdir(path)
os.makedirs('./Zn-Pc-allBr/'+str(xc),exist_ok=True)
os.chdir('./Zn-Pc-allBr/'+str(xc))
mollist=['Zn-Pc-allBr']
lat=[30.0,30.0,30.0]
bandscale=0
for mol in mollist:
    #if mol == 'Zn-Pc_1mer':
    #    bandscale=30
    #elif mol == 'Zn-Pc_2mer':
    #    bandscale=30
    os.makedirs(str(mol),exist_ok=True)
    os.chdir(str(mol))
    xyz = read(g.cifdir+'/'+str(mol)+'.xyz')
    xyz = sortmol(xyz)
    com = xyz.get_center_of_mass()
    shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
    #shift=[0.0,0.0,0.0]-com
    xyz.translate(shift)
    xyz.set_cell(lat)
    #view(xyz)
    for basis in ['DZP']:
        os.makedirs(str(basis),exist_ok=True)
        os.chdir(str(basis))
        for cutoff in [100.0]:
            os.makedirs('cutoff_'+str(cutoff),exist_ok=True)
            os.chdir('cutoff_'+str(cutoff))
            # no efield
            os.makedirs('e0',exist_ok=True)
            os.chdir('e0')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex+
            os.makedirs('ex+',exist_ok=True)
            os.chdir('ex+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex-
            os.makedirs('ex-',exist_ok=True)
            os.chdir('ex-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=-epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey+
            os.makedirs('ey+',exist_ok=True)
            os.chdir('ey+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey-
            os.makedirs('ey-',exist_ok=True)
            os.chdir('ey-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=-epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')            
            # ez+
            os.makedirs('ez+',exist_ok=True)
            os.chdir('ez+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=0.0,ez=epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ez-
            os.makedirs('ez-',exist_ok=True)
            os.chdir('ez-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,None,g.siesta_pot,bandscale=bandscale,ex=0.0,ey=0.0,ez=-epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')  
# 双極子モーメントと分極率テンソルの試計算終了
"""
            
"""
# Zn-フタロシアニンのパッキングモデル作成

#T_list = [300, 500, 700, 900,1100]
T_list = [1100]
nset = 1

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

os.makedirs('./Zn-Pc_train',exist_ok=True)
os.chdir('./Zn-Pc_train')

#os.makedirs('./Zn-Pc_md',exist_ok=True)
#os.chdir('./Zn-Pc_md')

os.makedirs('./Zn-Pc-allBr',exist_ok=True)
os.chdir('./Zn-Pc-allBr')

os.makedirs('./lammps',exist_ok=True)
os.chdir('./lammps')

xyz = read(g.cifdir+'/Zn-Pc-allBr.xyz')
#xyz = read(g.cifdir+'/Zn-Pc_1mer.xyz')

xyz = sort(xyz)

#xyz.write('mols.xyz')

mollist={
    'mols':64
}

x_box=100.0
y_box=100.0
z_box=100.0

for T0 in T_list:
    os.makedirs('./T_'+str(T0),exist_ok=True)
    os.chdir('./T_'+str(T0))
    for n in range(nset):
        os.makedirs('./set'+str(n),exist_ok=True)
        os.chdir('./set'+str(n))
        xyz.write('mols.xyz')
        mk_packmol_random(mollist,x_box,y_box,z_box)
        os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
        #
        cell = read('./system.xyz')
        cell.set_cell([x_box,y_box,z_box])
        cell = sortmol(cell)
        cell.write('system.cif')
        #
        #view(cell)
        #sys.exit()
        #
        unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
        #
        volume = unitcell.get_volume()
        masses = unitcell.get_masses()
        tot_mass = sum(masses)
        dens = tot_mass/volume
        print(dens)
        #
        volume1 = cell.get_volume()
        masses1 = cell.get_masses()
        tot_mass1 = sum(masses1)
        dens1 = tot_mass1/volume1
        print(dens1)
        #
        lat = (volume1/(dens/dens1))**(1/3)
        #lat = lat*(2.0)**(1/3)
        print(lat)
        #
        mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,200000,200000,T0,12345)
        os.chdir('../')
    os.chdir('../')

# Zn-フタロシアニンのパッキングモデル作成終了
"""


# Zn-フタロシアニンの学習データ作成

T_list = [300,500,700,900,1100]
nset = 1

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./Zn-Pc_train',exist_ok=True)
os.chdir('./Zn-Pc_train')

os.makedirs('./Zn-Pc-allBr',exist_ok=True)
os.chdir('./Zn-Pc-allBr')

os.makedirs('./siesta',exist_ok=True)
os.chdir('./siesta')

for xc in ['BLYP']:
    os.makedirs('./xc_'+str(xc),exist_ok=True)
    os.chdir('./xc_'+str(xc))
    for basis in ['DZP']:
        os.makedirs('./basis_'+str(basis),exist_ok=True)
        os.chdir('./basis_'+str(basis))
        for T0 in T_list:
            os.makedirs('./T_'+str(T0),exist_ok=True)
            os.chdir('./T_'+str(T0))
            for n in range(nset):
                os.makedirs('set'+str(n),exist_ok=True)
                os.chdir('set'+str(n))
                cell = read('/home/A23321P/work/myLAMMPS/dp_raman_test/Zn-Pc_train/Zn-Pc-allBr/lammps/T_'+str(T0)+'/set'+str(n)+'/result.data',format='lammps-data',sort_by_id=True,Z_of_type={1:35,2:6,3:7,4:30})
                #view(cell)
                #sys.exit()
                mk_siesta_input_npt(cell,xc,basis,100.0,None,T0,0.0,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
                os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
# Zn-フタロシアニンの学習データ作成終了


"""
# Zn-フタロシアニンの学習データ作成(リスタート)

T_list = [300,500,700,900,1100]
nset = 1

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./Zn-Pc_train',exist_ok=True)
os.chdir('./Zn-Pc_train')

os.makedirs('./Zn-Pc-allBr',exist_ok=True)
os.chdir('./Zn-Pc-allBr')

os.makedirs('./siesta',exist_ok=True)
os.chdir('./siesta')

for xc in ['KBM']:
    os.makedirs('./xc_'+str(xc),exist_ok=True)
    os.chdir('./xc_'+str(xc))
    for basis in ['DZP']:
        os.makedirs('./basis_'+str(basis),exist_ok=True)
        os.chdir('./basis_'+str(basis))
        for T0 in T_list:
            os.makedirs('./T_'+str(T0),exist_ok=True)
            os.chdir('./T_'+str(T0))
            for n in range(nset):
                os.makedirs('set'+str(n),exist_ok=True)
                os.chdir('set'+str(n))
                cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-allBr/datas/'+str(xc)+'/set3/siesta/xc_'+str(xc)+'/basis_'+str(basis)+'/T_'+str(T0)+'/set'+str(n)+'/siesta.STRUCT_OUT',format='struct_out')
                #view(cell)
                #sys.exit()
                mk_siesta_input_npt(cell,xc,basis,100.0,None,T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
                os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
# Zn-フタロシアニンの学習データ作成(リスタート)終了
"""

"""
# Zn-Pc Br置換モデル構造最適化
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./Zn-Pc_allBr',exist_ok=True)
os.chdir('./Zn-Pc_allBr')
os.makedirs('./1-mer',exist_ok=True)
os.chdir('./1-mer')
os.makedirs('./opt_Br',exist_ok=True)
os.chdir('./opt_Br')
xyz = read(g.cifdir+'/Zn-Pc_1mer.xyz')
symbols=[]
for atom in xyz:
    if atom.symbol == 'H':
        symbols.append('Br')
    else:
        symbols.append(atom.symbol)
xyz.set_chemical_symbols(symbols)
xyz = sortmol(xyz)
com = xyz.get_center_of_mass()
lat = [30.0,30.0,30.0]
shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
xyz.translate(shift)
xyz.set_cell(lat)
view(xyz)
mk_siesta_input_optimize(xyz,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
# Zn-Pc Br置換モデル構造最適化 終了
"""

"""
# Zn-Pc Br置換モデルWannier関数
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./Zn-Pc_allBr',exist_ok=True)
os.chdir('./Zn-Pc_allBr')
os.makedirs('./2-mer',exist_ok=True)
os.chdir('./2-mer')
os.makedirs('./wannier',exist_ok=True)
os.chdir('./wannier')
xyz = read(g.cifdir+'/Zn-Pc-allBr_2mer.xyz')
xyz = sortmol(xyz)
com = xyz.get_center_of_mass()
lat = [30.0,30.0,30.0]
shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
xyz.translate(shift)
xyz.set_cell(lat)
#view(xyz)
mk_siesta_input_scf_withEfield_wannier(xyz,'PBE','SZ',50.0,None,g.siesta_pot,bandscale=30,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
# Zn-Pc Br置換モデルWannier関数 終了
"""


"""
# Zn-フタロシアニン単分子構造最適化
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./ZincPhthalocyanine_test',exist_ok=True)
os.chdir('./ZincPhthalocyanine_test')
allBr = read(g.cifdir+'/ZincPhthalocyanine_allBr.xyz')
allBr_sorted = sort(allBr)
com = allBr_sorted.get_center_of_mass()
lat = [30,30,30]
shift = [15,15,15]-com
allBr_sorted.translate(shift)
allBr_sorted.set_cell(lat)
mk_siesta_input_optimize(allBr_sorted,'VDW','DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
#view(allBr)
# Zn-フタロシアニン単分子構造最適化終了
"""

"""
# H2O 単分子テスト
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./h2o_monomer',exist_ok=True)
os.chdir('./h2o_monomer')
h2o_mol = molecule('H2O')
com = h2o_mol.get_center_of_mass()
shift = [5.0,5.0,5.0] - com
h2o_mol.translate(shift)
h2o_mol.set_cell([10.0,10.0,10.0])
view(h2o_mol)
#mk_siesta_input_optimize(h2o_mol,'VDW','SZ',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
mk_siesta_input_nvt(h2o_mol,'VDW','SZ',50.0,[1,1,1],300.0,10,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
# H2O 単分子テスト終了
"""