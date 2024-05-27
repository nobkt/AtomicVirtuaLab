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
from AtomicVirtuaLab.tools import add_displacement, deform_cell
from AtomicVirtuaLab.io import cell2atomlist
from AtomicVirtuaLab.gpaw import mk_gpaw_lcao_input_optimize,mk_gpaw_pw_input_optimize,mk_gpaw_pw_input_npt,mk_gpaw_lcao_input_npt
from ase.data import atomic_numbers
import shutil
import os
import sys
import random

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"



#
# 分極率学習データ変換
#
#
# Wannier XYZ読み込み
#
os.makedirs('dp_pol_test',exist_ok=True)
os.chdir('dp_pol_test')

dir_ = '/home/A23321P/work/mySiesta/dipole_and_polarizability/wannier/Zn-Pc_train_3'
#dir_ = '/home/A23321P/work/mySiesta/dipole_and_polarizability/wannier/H-Zn-Pc_train_r4'
l_T = [300,500,700,900,1100]
l_set = [1]
l_e = ['ex-','ex+','ey-','ey+','ez-','ez+']
for T in l_T:
    os.chdir(dir_+'/T_'+str(T))
    for iset in l_set:
        os.chdir('set'+str(iset))
        os.chdir('traj1999')
        for e in l_e:
            os.chdir(e)
            #xyz = read('siesta_centres.xyz')
            xyz = read('siesta.xyz')
            symbols = xyz.get_chemical_symbols()
            #print(symbols)
            #tmp=[]
            #for symbol in symbols:
            #    if symbol == 'X':
            #        tmp.append('O')
            #    else:
            #        tmp.append(symbol)
            #xyz.set_chemical_symbols(tmp)
            lat = xyz.get_cell()
            #view(xyz)
            f = open('siesta.win','r')
            lines = f.readlines()
            f.close()
            i = 0
            for line in lines:
                if 'begin unit_cell_cart' in line:
                    lat[0][0] = float(lines[i+2].split()[0])
                    lat[0][1] = float(lines[i+2].split()[1])
                    lat[0][2] = float(lines[i+2].split()[2])
                    lat[1][0] = float(lines[i+3].split()[0])
                    lat[1][1] = float(lines[i+3].split()[1])
                    lat[1][2] = float(lines[i+3].split()[2])
                    lat[2][0] = float(lines[i+4].split()[0])
                    lat[2][1] = float(lines[i+4].split()[1])
                    lat[2][2] = float(lines[i+4].split()[2])                     
                i = i + 1
            xyz.set_cell(lat)
            xyz.set_pbc([1,1,1])
            xyz.write('POSCAR')
            cell = sortmol(xyz,sort_atom=True)
            #view(cell)
            cell.write('test.xyz')
            #cell.write('POSCAR')
            sys.exit()
           
            
        
        
    

"""
# MD 比率
subelm = 'H'
sublists = []
sublists.append({
    'slist' : {'Cl':0,'Br':16},
    'nmol' : 12    
})

sublists.append({
    'slist' : {'Cl':1,'Br':15},
    'nmol' : 14    
})

sublists.append({
    'slist' : {'Cl':0,'Br':15},
    'nmol' : 3    
})

sublists.append({
    'slist' : {'Cl':2,'Br':14},
    'nmol' : 10    
})

sublists.append({
    'slist' : {'Cl':1,'Br':14},
    'nmol' : 4    
})

sublists.append({
    'slist' : {'Cl':3,'Br':13},
    'nmol' : 6    
})

sublists.append({
    'slist' : {'Cl':2,'Br':13},
    'nmol' : 3    
})

sublists.append({
    'slist' : {'Cl':4,'Br':12},
    'nmol' : 3    
})

sublists.append({
    'slist' : {'Cl':3,'Br':12},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':2,'Br':12},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':5,'Br':11},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':4,'Br':11},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':6,'Br':10},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':5,'Br':10},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':7,'Br':9},
    'nmol' : 1    
})

nmols = 64
"""

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

"""
os.makedirs('A110',exist_ok=True)
os.chdir('./A110')

xyz = read(g.cifdir+'/Zn-Pc-allH.xyz')
xyz = sort(xyz)

mollist={
    'mols':nmols
}

x_box=100.0
y_box=100.0
z_box=100.0

os.makedirs('./tmp',exist_ok=True)
os.chdir('./tmp')
xyz.write('mols.xyz')
mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('./system.xyz')
cell.set_cell([x_box,y_box,z_box])
cell = sortmol(cell)
symbols=[]
nlist=0
nlist_tmp1=1
nlist_tmp2=0
for atom in cell:
    symbols.append(atom.symbol)
for imol in range(nmols):
    ltmp=[]
    aid=0
    for symbol in symbols:
        if cell.arrays['mol-id'][aid] == imol+1 and symbol == subelm:
            ltmp.append(aid)
        aid=aid+1
    #print(ltmp)
    nsub=0
    subnum={}
    sublist = sublists[nlist]
    if nlist_tmp1 > nlist_tmp2 + sublist['nmol']:
        nlist_tmp2 = nlist_tmp2 + sublist['nmol']
        nlist = nlist + 1
        sublist = sublists[nlist]
    for elm0 in sublist['slist']:
        #n0 = random.randint(0, sublist[elm0])
        n0 = sublist['slist'][elm0]
        subnum[elm0] = n0
        nsub = nsub + n0
    subs = random.sample(ltmp, nsub)
    nn = 0
    for elm0 in subnum:
        n0 = subnum[elm0]
        for i in range(nn,nn+n0):
            cell[subs[i]].symbol = elm0
        nn = nn + n0
    nlist_tmp1 = nlist_tmp1 + 1

#cell = sortmol(cell)
cell.write('system.cif')
symbols = cell2atomlist(cell)
#
unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
#
volume = unitcell.get_volume()
masses = unitcell.get_masses()
tot_mass = sum(masses)
dens = tot_mass/volume
#
volume1 = cell.get_volume()
masses1 = cell.get_masses()
tot_mass1 = sum(masses1)
dens1 = tot_mass1/volume1
#
lat = (volume1/(dens/dens1))**(1/3)
#
mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,20000,2000,300.0,12345)
os.system('mpirun -np 20 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
Z_of_type={}
i = 0
for symbol in symbols:
    Z_of_type[i+1] = atomic_numbers[symbol]
    i=i+1
cell = read('result.data',format='lammps-data',sort_by_id=True,Z_of_type=Z_of_type)
#view(cell)
os.chdir('../')
os.makedirs('npt',exist_ok=True)
os.chdir('npt')
mk_npt_input_deepmd(cell,0.0005,2000,2000,2000000,300,0.0,12345,mol=True)
#os.chdir('opt')
#mk_siesta_input_optimize(cell,'PBE','DZP',50.0,[1,1,1],50,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized')
##os.chdir('../')
#os.makedirs('npt',exist_ok=True)
#os.chdir('npt')
#mk_siesta_input_npt(cell,'VDW','DZP',100,[1,1,1],T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
#os.chdir('../')
#mk_gpaw_pw_input_optimize(cell,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#mk_gpaw_pw_input_npt(cell,T0,0.0,2000,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=400,restart=True)
#mk_gpaw_lcao_input_optimize(cell,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#sys.exit()
os.chdir('../')

"""

"""
# MD 比率
subelm = 'H'
sublists = []
sublists.append({
    'slist' : {'Cl':0,'Br':16},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':1,'Br':15},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':0,'Br':15},
    'nmol' : 4    
})

sublists.append({
    'slist' : {'Cl':2,'Br':14},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':1,'Br':14},
    'nmol' : 4    
})

sublists.append({
    'slist' : {'Cl':0,'Br':14},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':3,'Br':13},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':2,'Br':13},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':1,'Br':13},
    'nmol' : 6    
})

sublists.append({
    'slist' : {'Cl':0,'Br':13},
    'nmol' : 4    
})

sublists.append({
    'slist' : {'Cl':3,'Br':12},
    'nmol' : 4    
})

sublists.append({
    'slist' : {'Cl':2,'Br':12},
    'nmol' : 6    
})

sublists.append({
    'slist' : {'Cl':1,'Br':12},
    'nmol' : 5    
})

sublists.append({
    'slist' : {'Cl':0,'Br':12},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':5,'Br':11},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':4,'Br':11},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':3,'Br':11},
    'nmol' : 3    
})

sublists.append({
    'slist' : {'Cl':2,'Br':11},
    'nmol' : 5    
})

sublists.append({
    'slist' : {'Cl':6,'Br':10},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':5,'Br':10},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':4,'Br':10},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':7,'Br':9},
    'nmol' : 2    
})

sublists.append({
    'slist' : {'Cl':6,'Br':9},
    'nmol' : 1    
})

sublists.append({
    'slist' : {'Cl':8,'Br':8},
    'nmol' : 3    
})

nmols = 64

os.makedirs('C100',exist_ok=True)
os.chdir('./C100')

xyz = read(g.cifdir+'/Zn-Pc-allH.xyz')
xyz = sort(xyz)

mollist={
    'mols':nmols
}

x_box=100.0
y_box=100.0
z_box=100.0

os.makedirs('./tmp',exist_ok=True)
os.chdir('./tmp')
xyz.write('mols.xyz')
mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('./system.xyz')
cell.set_cell([x_box,y_box,z_box])
cell = sortmol(cell)
symbols=[]
nlist=0
nlist_tmp1=1
nlist_tmp2=0
for atom in cell:
    symbols.append(atom.symbol)
for imol in range(nmols):
    ltmp=[]
    aid=0
    for symbol in symbols:
        if cell.arrays['mol-id'][aid] == imol+1 and symbol == subelm:
            ltmp.append(aid)
        aid=aid+1
    #print(ltmp)
    nsub=0
    subnum={}
    sublist = sublists[nlist]
    if nlist_tmp1 > nlist_tmp2 + sublist['nmol']:
        nlist_tmp2 = nlist_tmp2 + sublist['nmol']
        nlist = nlist + 1
        sublist = sublists[nlist]
    for elm0 in sublist['slist']:
        #n0 = random.randint(0, sublist[elm0])
        n0 = sublist['slist'][elm0]
        subnum[elm0] = n0
        nsub = nsub + n0
    subs = random.sample(ltmp, nsub)
    nn = 0
    for elm0 in subnum:
        n0 = subnum[elm0]
        for i in range(nn,nn+n0):
            cell[subs[i]].symbol = elm0
        nn = nn + n0
    nlist_tmp1 = nlist_tmp1 + 1

#cell = sortmol(cell)
cell.write('system.cif')
symbols = cell2atomlist(cell)
#
unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
#
volume = unitcell.get_volume()
masses = unitcell.get_masses()
tot_mass = sum(masses)
dens = tot_mass/volume
#
volume1 = cell.get_volume()
masses1 = cell.get_masses()
tot_mass1 = sum(masses1)
dens1 = tot_mass1/volume1
#
lat = (volume1/(dens/dens1))**(1/3)
#
mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,20000,2000,300.0,12345)
os.system('mpirun -np 20 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
Z_of_type={}
i = 0
for symbol in symbols:
    Z_of_type[i+1] = atomic_numbers[symbol]
    i=i+1
cell = read('result.data',format='lammps-data',sort_by_id=True,Z_of_type=Z_of_type)
#view(cell)
os.chdir('../')
os.makedirs('npt',exist_ok=True)
os.chdir('npt')
mk_npt_input_deepmd(cell,0.0005,2000,2000,2000000,300,0.0,12345,mol=True)
#os.chdir('opt')
#mk_siesta_input_optimize(cell,'PBE','DZP',50.0,[1,1,1],50,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized')
##os.chdir('../')
#os.makedirs('npt',exist_ok=True)
#os.chdir('npt')
#mk_siesta_input_npt(cell,'VDW','DZP',100,[1,1,1],T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
#os.chdir('../')
#mk_gpaw_pw_input_optimize(cell,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#mk_gpaw_pw_input_npt(cell,T0,0.0,2000,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=400,restart=True)
#mk_gpaw_lcao_input_optimize(cell,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#sys.exit()
os.chdir('../')
"""

"""
# MD 
subelm = 'Cl'
sublist = {'H':0,'Br':15}
nmols = 64

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

os.makedirs('Br15Cl1',exist_ok=True)
os.chdir('./Br15Cl1')

xyz = read(g.cifdir+'/Zn-Pc-allH.xyz')
xyz = sort(xyz)

mollist={
    'mols':nmols
}

x_box=100.0
y_box=100.0
z_box=100.0

os.makedirs('./tmp',exist_ok=True)
os.chdir('./tmp')
xyz.write('mols.xyz')
mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('./system.xyz')
cell.set_cell([x_box,y_box,z_box])
cell = sortmol(cell)
symbols=[]
for atom in cell:
    symbols.append(atom.symbol)
for imol in range(nmols):
    ltmp=[]
    aid=0
    for symbol in symbols:
        if cell.arrays['mol-id'][aid] == imol+1 and symbol == subelm:
            ltmp.append(aid)
        aid=aid+1
    #print(ltmp)
    nsub=0
    subnum={}
    for elm0 in sublist:
        #n0 = random.randint(0, sublist[elm0])
        n0 = sublist[elm0]
        subnum[elm0] = n0
        nsub = nsub + n0
    subs = random.sample(ltmp, nsub)
    nn = 0
    for elm0 in subnum:
        n0 = subnum[elm0]
        for i in range(nn,nn+n0):
            cell[subs[i]].symbol = elm0
        nn = nn + n0
#cell = sortmol(cell)
cell.write('system.cif')
symbols = cell2atomlist(cell)
#
unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
#
volume = unitcell.get_volume()
masses = unitcell.get_masses()
tot_mass = sum(masses)
dens = tot_mass/volume
#
volume1 = cell.get_volume()
masses1 = cell.get_masses()
tot_mass1 = sum(masses1)
dens1 = tot_mass1/volume1
#
lat = (volume1/(dens/dens1))**(1/3)
#
mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,20000,2000,300.0,12345)
os.system('mpirun -np 20 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
Z_of_type={}
i = 0
for symbol in symbols:
    Z_of_type[i+1] = atomic_numbers[symbol]
    i=i+1
cell = read('result.data',format='lammps-data',sort_by_id=True,Z_of_type=Z_of_type)
#view(cell)
os.chdir('../')
os.makedirs('npt',exist_ok=True)
os.chdir('npt')
mk_npt_input_deepmd(cell,0.0005,2000,2000,2000000,300,0.0,12345,mol=True)
#os.chdir('opt')
#mk_siesta_input_optimize(cell,'PBE','DZP',50.0,[1,1,1],50,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized')
##os.chdir('../')
#os.makedirs('npt',exist_ok=True)
#os.chdir('npt')
#mk_siesta_input_npt(cell,'VDW','DZP',100,[1,1,1],T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
#os.chdir('../')
#mk_gpaw_pw_input_optimize(cell,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#mk_gpaw_pw_input_npt(cell,T0,0.0,2000,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=400,restart=True)
#mk_gpaw_lcao_input_optimize(cell,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
#sys.exit()
os.chdir('../')
"""

"""
# Zn-フタロシアニンの学習データ作成 Siesta

T_list = [300, 500, 700, 900, 1100]
#T_list = [1100]
nset = 10
nmols = 6
subelm = 'X'
#sublist = {'Cl':4,'H':4}
sublist = {}

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

os.makedirs('./Zn-Pc_train',exist_ok=True)
os.chdir('./Zn-Pc_train')

xyz = read(g.cifdir+'/Zn-Pc-allCl.xyz')
xyz = sort(xyz)

mollist={
    'mols':nmols
}

x_box=100.0
y_box=100.0
z_box=100.0

for n in range(nset):
    for T0 in T_list:
        os.makedirs('./T_'+str(T0),exist_ok=True)
        os.chdir('./T_'+str(T0))
        os.makedirs('./set'+str(n),exist_ok=True)
        os.chdir('./set'+str(n))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        xyz.write('mols.xyz')
        mk_packmol_random(mollist,x_box,y_box,z_box)
        os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
        #
        cell = read('./system.xyz')
        cell.set_cell([x_box,y_box,z_box])
        cell = sortmol(cell)
        symbols=[]
        for atom in cell:
            symbols.append(atom.symbol)
        for imol in range(nmols):
            ltmp=[]
            aid=0
            for symbol in symbols:
                if cell.arrays['mol-id'][aid] == imol+1 and symbol == subelm:
                    ltmp.append(aid)
                aid=aid+1
            #print(ltmp)
            nsub=0
            subnum={}
            for elm0 in sublist:
                n0 = random.randint(0, sublist[elm0])
                subnum[elm0] = n0
                nsub = nsub + n0
            subs = random.sample(ltmp, nsub)
            nn = 0
            for elm0 in subnum:
                n0 = subnum[elm0]
                for i in range(nn,nn+n0):
                    cell[subs[i]].symbol = elm0
                nn = nn + n0
        #cell = sortmol(cell)
        cell.write('system.cif')
        symbols = cell2atomlist(cell)
        #
        unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
        #
        volume = unitcell.get_volume()
        masses = unitcell.get_masses()
        tot_mass = sum(masses)
        dens = tot_mass/volume
        #
        volume1 = cell.get_volume()
        masses1 = cell.get_masses()
        tot_mass1 = sum(masses1)
        dens1 = tot_mass1/volume1
        #
        lat = (volume1/(dens/dens1))**(1/3)
        #
        mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,20000,2000,T0,12345)
        os.system('mpirun -np 2 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
        Z_of_type={}
        i = 0
        for symbol in symbols:
            Z_of_type[i+1] = atomic_numbers[symbol]
            i=i+1
        cell = read('result.data',format='lammps-data',sort_by_id=True,Z_of_type=Z_of_type)
        #view(cell)
        os.chdir('../')
        os.makedirs('opt',exist_ok=True)
        os.chdir('opt')
        mk_siesta_input_optimize(cell,'PBE','DZP',50.0,[1,1,1],50,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized')
        os.chdir('../')
        os.makedirs('npt',exist_ok=True)
        os.chdir('npt')
        mk_siesta_input_npt(cell,'VDW','DZP',100,[1,1,1],T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
        os.chdir('../')
        #mk_gpaw_pw_input_optimize(cell,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
        #mk_gpaw_pw_input_npt(cell,T0,0.0,2000,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=400,restart=True)
        #mk_gpaw_lcao_input_optimize(cell,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
        #sys.exit()
        os.chdir('../')
        os.chdir('../')

# Zn-フタロシアニンの学習データ作成 Siesta
"""


"""
# Zn-フタロシアニンの学習データ作成 GPAW

T_list = [300, 500, 700, 900, 1100]
#T_list = [1100]
nset = 10
nmols = 6
subelm = 'Br'
sublist = {'Cl':4,'H':4}

os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')

os.makedirs('./Zn-Pc_train',exist_ok=True)
os.chdir('./Zn-Pc_train')

xyz = read(g.cifdir+'/Zn-Pc-allBr.xyz')
xyz = sort(xyz)

mollist={
    'mols':nmols
}

x_box=100.0
y_box=100.0
z_box=100.0

for n in range(nset):
    for T0 in T_list:
        os.makedirs('./T_'+str(T0),exist_ok=True)
        os.chdir('./T_'+str(T0))
        os.makedirs('./set'+str(n),exist_ok=True)
        os.chdir('./set'+str(n))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        xyz.write('mols.xyz')
        mk_packmol_random(mollist,x_box,y_box,z_box)
        os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
        #
        cell = read('./system.xyz')
        cell.set_cell([x_box,y_box,z_box])
        cell = sortmol(cell)
        symbols=[]
        for atom in cell:
            symbols.append(atom.symbol)
        for imol in range(nmols):
            ltmp=[]
            aid=0
            for symbol in symbols:
                if cell.arrays['mol-id'][aid] == imol+1 and symbol == subelm:
                    ltmp.append(aid)
                aid=aid+1
            #print(ltmp)
            nsub=0
            subnum={}
            for elm0 in sublist:
                n0 = random.randint(0, sublist[elm0])
                subnum[elm0] = n0
                nsub = nsub + n0
            subs = random.sample(ltmp, nsub)
            nn = 0
            for elm0 in subnum:
                n0 = subnum[elm0]
                for i in range(nn,nn+n0):
                    cell[subs[i]].symbol = elm0
                nn = nn + n0
        cell.write('system.cif')
        symbols = cell2atomlist(cell)
        #
        unitcell = read(g.cifdir+'/ZincPhthalocyanine.cif')
        #
        volume = unitcell.get_volume()
        masses = unitcell.get_masses()
        tot_mass = sum(masses)
        dens = tot_mass/volume
        #
        volume1 = cell.get_volume()
        masses1 = cell.get_masses()
        tot_mass1 = sum(masses1)
        dens1 = tot_mass1/volume1
        #
        lat = (volume1/(dens/dens1))**(1/3)
        #
        mk_nvt_input_uff_rigid_scale(cell,0.0005,1000,1000,lat,20000,2000,T0,12345)
        os.system('mpirun -np 1 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
        Z_of_type={}
        i = 0
        for symbol in symbols:
            Z_of_type[i+1] = atomic_numbers[symbol]
            i=i+1
        cell = read('result.data',format='lammps-data',sort_by_id=True,Z_of_type=Z_of_type)
        #view(cell)
        os.chdir('../')
        #mk_gpaw_pw_input_optimize(cell,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
        #mk_gpaw_pw_input_npt(cell,T0,0.0,2000,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=400,restart=True)
        mk_gpaw_lcao_input_optimize(cell,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=True)
        mk_gpaw_lcao_input_npt(cell,T0,0.0,500,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=True,berendsen=True,berandsen_nstep=100,restart=True)
        #sys.exit()
        os.chdir('../')
        os.chdir('../')

# Zn-フタロシアニンの学習データ作成 GPAW
"""


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
dpdata=False
if dpdata:
    path = os.getcwd()
    datadir_ = '/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc_tmp'
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
# Wannier関数学習データ
au2v=51.4220632
epsilon=0.001*au2v
xc='VDW'
cutoff=100.0
basis='DZP'
os.makedirs('./dp_raman_test',exist_ok=True)
os.chdir('./dp_raman_test')
os.makedirs('./dipole_and_polarizability',exist_ok=True)
os.chdir('./dipole_and_polarizability')

path = os.getcwd()
os.chdir(path)

# Br-Cl-H_Zn-Pc
os.chdir(path)
mollist=['Zn-Pc_train_3']
Tlist=[300,500,700,900,1100]
nset=10
bandscale=30
for mol in mollist:
    os.makedirs(mol,exist_ok=True)
    os.chdir('./'+mol)
    for T0 in Tlist:
        os.makedirs('T_'+str(T0),exist_ok=True)
        os.chdir('T_'+str(T0))
        for iset in range(nset):
            os.makedirs('set'+str(iset),exist_ok=True)
            os.chdir('set'+str(iset))
            cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc_tmp/'+mol+'/T_'+str(T0)+'/set'+str(iset)+'/npt/siesta.AXSF',format='xsf',index=':')
            #view(cell)
            #sys.exit()
            itraj = 0
            for xyz in cell:
                if itraj < 1499:
                    itraj = itraj + 1
                    continue
                os.makedirs('traj'+str(itraj),exist_ok=True)
                os.chdir('traj'+str(itraj))
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
                #sys.exit()
                itraj = itraj + 1
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

# M_Zn-Pc
os.chdir(path)
mollist=['H-Zn-Pc_train_r4','Cl-Zn-Pc_train_r4','Br-Zn-Pc_train_r4']
Tlist=[300,500,700,900,1100]
nset=1
bandscale=30
for mol in mollist:
    os.makedirs(mol,exist_ok=True)
    os.chdir('./'+mol)
    for T0 in Tlist:
        os.makedirs('T_'+str(T0),exist_ok=True)
        os.chdir('T_'+str(T0))
        for iset in range(nset):
            os.makedirs('set'+str(iset),exist_ok=True)
            os.chdir('set'+str(iset))
            cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc_tmp/'+mol+'/T_'+str(T0)+'/set'+str(iset)+'/npt/siesta.AXSF',format='xsf',index=':')
            #view(cell)
            #sys.exit()
            itraj = 0
            for xyz in cell:
                if itraj < 1499:
                    itraj = itraj + 1
                    continue
                os.makedirs('traj'+str(itraj),exist_ok=True)
                os.chdir('traj'+str(itraj))
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
                #sys.exit()
                itraj = itraj + 1
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

# Wannier関数学習データ 終了
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

for xc in ['PBE']:
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
                mk_siesta_input_npt(cell,xc,basis,200.0,None,T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
                os.chdir('../')
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
# Zn-フタロシアニンの学習データ作成終了
"""

"""
# Zn-フタロシアニンの学習データ作成(リスタート)

T_list = [300,500,700,900,1100]
nset = [1]

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
            for n in nset:
                os.makedirs('set'+str(n),exist_ok=True)
                os.chdir('set'+str(n))
                cell = read('/home/A23321P/work/mySiesta/dipole_and_polarizability/Zn-Pc-allBr/datas/'+str(xc)+'/xc_'+str(xc)+'/basis_'+str(basis)+'/T_'+str(T0)+'/set0/siesta.STRUCT_OUT',format='struct_out')
                cell = add_displacement(cell,0.05)
                cell = deform_cell(cell,0.05,0.5)
                #view(cell)
                #sys.exit()
                mk_siesta_input_npt(cell,xc,basis,200.0,None,T0,0.0,500,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
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
os.makedirs('opt',exist_ok=True)
os.chdir('opt')
mk_siesta_input_optimize(h2o_mol,'VDW','SZ',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
os.chdir('../')
os.makedirs('nvt',exist_ok=True)
os.chdir('nvt')
mk_siesta_input_nvt(h2o_mol,'VDW','SZ',50.0,[1,1,1],300.0,10,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized')
os.chdir('../')
# H2O 単分子テスト終了
"""