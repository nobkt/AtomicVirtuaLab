from ase.build import molecule
from ase.io import read, write
from ase.visualize import view
from ase.build import sort, make_supercell
from AtomicVirtuaLab.espresso import mk_qe_input_relax
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize,mk_siesta_input_nvt,mk_siesta_input_scf_withEfield,mk_siesta_input_scf_withEfield_wannier,mk_siesta_input_cellopt
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_nvt_input_uff_rigid
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'


# Zn-Pc アモルファスモデル読み込み
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('Zn-Pc-amo',exist_ok=True)
os.chdir('Zn-Pc-amo')

cell = read(g.cifdir+'/Zn-Pc_10.data',format='lammps-data',sort_by_id=True,Z_of_type={1:6,2:1,3:7,4:30})

view(cell)

#mk_siesta_input_optimize(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

mk_siesta_input_cellopt(cell,'PBE','SZ',50.0,None,2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')

# Zn-Pc アモルファスモデル読み込み 終了

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

# Zn-Pc 結晶モデル作成
"""


"""
# 双極子モーメントと分極率テンソルの試計算
au2v=51.4220632
epsilon=0.001*au2v
xc='PBE'
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./dipole_and_polarizability',exist_ok=True)
os.chdir('./dipole_and_polarizability')

path = os.getcwd()
os.chdir(path)

# 水クラスター
os.makedirs('./water_clusters/'+str(xc),exist_ok=True)
os.chdir('./water_clusters/'+str(xc))
mollist=['H2O_1mer','H2O_5mer','H2O_10mer','H2O_20mer']
lat=[30.0,30.0,30.0]
for mol in mollist:
    if mol == 'H2O_1mer':
        bandscale = 0
    else:
        bandscale = 0
    os.makedirs(str(mol),exist_ok=True)
    os.chdir(str(mol))
    xyz = read(g.cifdir+'/'+str(mol)+'.xyz')
    xyz = sortmol(xyz)
    com = xyz.get_center_of_mass()
    shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
    xyz.translate(shift)
    xyz.set_cell(lat)
    #view(xyz)
    for basis in ['SZ']:
        os.makedirs(str(basis),exist_ok=True)
        os.chdir(str(basis))
        for cutoff in [50.0]:
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
# 水クラスター終了

# Zn-フタロシアニン
os.chdir(path)
os.makedirs('./Zn-Pc/'+str(xc),exist_ok=True)
os.chdir('./Zn-Pc/'+str(xc))
mollist=['Zn-Pc_1mer','Zn-Pc_2mer']
lat=[30.0,30.0,30.0]
for mol in mollist:
    if mol == 'Zn-Pc_1mer':
        bandscale=30
    elif mol == 'Zn-Pc_2mer':
        bandscale=30
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
    for basis in ['SZ']:
        os.makedirs(str(basis),exist_ok=True)
        os.chdir(str(basis))
        for cutoff in [50.0]:
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
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./ZincPhthalocyanine_2mol_test',exist_ok=True)
os.chdir('./ZincPhthalocyanine_2mol_test')
allBr = read(g.cifdir+'/Zn-Pc_1mer.xyz')
allBr_sorted = sort(allBr)

allBr_sorted.write('./test.xyz')

mollist={
    'test':15
}

x_box=50.0
y_box=50.0
z_box=50.0

mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('./system.xyz')
cell.set_cell([x_box,y_box,z_box])
cell = sortmol(cell)
cell.write('test.cif')

#mk_siesta_input_scf_withEfield_wannier(cell,'PBE','SZ',50.0,None,g.siesta_pot,bandscale=30,ex=0.0,ey=0.0,ez=-0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
#mk_qe_input_relax(cell,'pbe','us',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=36,ecutrho=400,mixing_beta=0.2,kpts=None,ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
mk_nvt_input_uff_rigid(cell,0.0005,100,100,20000,300.0,12345)
view(cell)

unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')

volume = unitcell.get_volume()
masses = unitcell.get_masses()
tot_mass = sum(masses)
dens = tot_mass/volume
print(dens)

volume1 = cell.get_volume()
masses1 = cell.get_masses()
tot_mass1 = sum(masses1)
dens1 = tot_mass1/volume1
print(dens1)

scale = (volume1/(dens/dens1))**(1/3)
print(scale)

# Zn-フタロシアニンのパッキングモデル作成終了
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