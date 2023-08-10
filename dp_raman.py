from ase.build import molecule
from ase.io import read, write
from ase.visualize import view
from ase.build import sort
from AtomicVirtuaLab.espresso import mk_qe_input_relax
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize,mk_siesta_input_nvt,mk_siesta_input_scf_withEfield,mk_siesta_input_scf_withEfield_wannier
from AtomicVirtuaLab.packmol import mk_packmol_random
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'


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
os.makedirs('./water_clusters',exist_ok=True)
os.chdir('./water_clusters')
mollist=['H2O_1mer','H2O_5mer','H2O_10mer','H2O_20mer']
lat=[30.0,30.0,30.0]
for mol in mollist:
    os.makedirs(str(mol),exist_ok=True)
    os.chdir(str(mol))
    xyz = read(g.cifdir+'/'+str(mol)+'.xyz')
    xyz = sortmol(xyz)
    com = xyz.get_center_of_mass()
    shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
    xyz.translate(shift)
    xyz.set_cell(lat)
    #view(xyz)
    for basis in ['DZ','DZP']:
        os.makedirs(str(basis),exist_ok=True)
        os.chdir(str(basis))
        for cutoff in [81.0]:
            os.makedirs('cutoff_'+str(cutoff),exist_ok=True)
            os.chdir('cutoff_'+str(cutoff))
            # no efield
            os.makedirs('e0',exist_ok=True)
            os.chdir('e0')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex+
            os.makedirs('ex+',exist_ok=True)
            os.chdir('ex+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex-
            os.makedirs('ex-',exist_ok=True)
            os.chdir('ex-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=-epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey+
            os.makedirs('ey+',exist_ok=True)
            os.chdir('ey+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey-
            os.makedirs('ey-',exist_ok=True)
            os.chdir('ey-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=-epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')            
            # ez+
            os.makedirs('ez+',exist_ok=True)
            os.chdir('ez+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ez-
            os.makedirs('ez-',exist_ok=True)
            os.chdir('ez-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=-epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')                   
# 水クラスター終了

# Zn-フタロシアニン
os.chdir(path)
os.makedirs('./Zn-Pc',exist_ok=True)
os.chdir('./Zn-Pc')
mollist=['Zn-Pc_1mer','Zn-Pc_2mer']
lat=[30.0,30.0,30.0]
for mol in mollist:
    os.makedirs(str(mol),exist_ok=True)
    os.chdir(str(mol))
    xyz = read(g.cifdir+'/'+str(mol)+'.xyz')
    xyz = sortmol(xyz)
    com = xyz.get_center_of_mass()
    shift=[lat[0]/2.0,lat[1]/2.0,lat[2]/2.0]-com
    xyz.translate(shift)
    xyz.set_cell(lat)
    #view(xyz)
    for basis in ['DZ','DZP']:
        os.makedirs(str(basis),exist_ok=True)
        os.chdir(str(basis))
        for cutoff in [81.0]:
            os.makedirs('cutoff_'+str(cutoff),exist_ok=True)
            os.chdir('cutoff_'+str(cutoff))
            # no efield
            os.makedirs('e0',exist_ok=True)
            os.chdir('e0')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex+
            os.makedirs('ex+',exist_ok=True)
            os.chdir('ex+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ex-
            os.makedirs('ex-',exist_ok=True)
            os.chdir('ex-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=-epsilon,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey+
            os.makedirs('ey+',exist_ok=True)
            os.chdir('ey+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ey-
            os.makedirs('ey-',exist_ok=True)
            os.chdir('ey-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=-epsilon,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')            
            # ez+
            os.makedirs('ez+',exist_ok=True)
            os.chdir('ez+')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            # ez-
            os.makedirs('ez-',exist_ok=True)
            os.chdir('ez-')            
            mk_siesta_input_scf_withEfield_wannier(xyz,xc,basis,cutoff,[1,1,1],g.siesta_pot,ex=0.0,ey=0.0,ez=-epsilon,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
            os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')  
# 双極子モーメントと分極率テンソルの試計算終了

            
"""
# Zn-フタロシアニンのパッキングモデル作成
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./ZincPhthalocyanine_2mol_test',exist_ok=True)
os.chdir('./ZincPhthalocyanine_2mol_test')
allBr = read(g.cifdir+'/ZincPhthalocyanine_allBr.xyz')
allBr_sorted = sort(allBr)

allBr_sorted.write('./test.xyz')

mollist={
    'test':2
}

x_box=30.0
y_box=30.0
z_box=30.0

mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('./system.xyz')
cell.set_cell([x_box,y_box,z_box])
cell.write('test.cif')

mk_siesta_input_optimize(cell,'VDW','DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
#mk_qe_input_relax(cell,'pbe','us',level='high',estep=1000,nstep=1000,nosym=True,ecutwfc=36,ecutrho=400,mixing_beta=0.2,kpts=None,ecut='manual',tstress=False,options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=False)
view(cell)

# Zn-フタロシアニンのパッキングモデル作成終了
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