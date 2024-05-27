from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.build import make_supercell
from ase.constraints import FixAtoms
from AtomicVirtuaLab.build import slabgen,sortmol
from ase.io import read
from AtomicVirtuaLab.espresso import mk_qe_input_relax
from AtomicVirtuaLab.lammps import mk_nvt_input_deepmd_friction, mk_opt_input_deepmd, mk_mimize_input_uff_adsorp_onMoS2, mk_mimize_input_dp_adsorp_onMoS2
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize, mk_siesta_input_scf
import shutil
import os
import sys
import random
import pandas as pd

g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"
g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'

os.makedirs('./MoS2_friction',exist_ok=True)
os.chdir('./MoS2_friction')



"""
cell1 = rd_cif(g.cifdir+'/MoS2_mp1434.cif')
cell2 = rd_cif(g.cifdir+'/MoS2_mp2815.cif')


# 6-layers
nx=9
ny=5
nz1=2
nz2=3
nlayer1=6
nlayer2=6

# 6-layers xy expand
#nx=18
#ny=10
#nz1=2
#nz2=3
#nlayer1=6
#nlayer2=6

# 20-layers
#nx=9
#ny=5
#nz1=7
#nz2=10
#nlayer1=21
#nlayer2=20



cell1 = slabgen(cell1,0,0,1,1,1,nz1,10.0,10.0)
os.system('atomsk slab3.cif -orthogonal-cell MoS2_mp1434_slab_3_ortho.cfg')
#view(cell1)
cell1 = read('MoS2_mp1434_slab_3_ortho.cfg')

cell2 = slabgen(cell2,0,0,1,1,1,nz2,10.0,10.0)
#view(cell2)
os.system('atomsk slab1.cif -orthogonal-cell MoS2_mp2815_slab_1_ortho.cfg')
cell2 = read('MoS2_mp2815_slab_1_ortho.cfg')

#cell1 = make_supercell(cell1,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
#cell2 = make_supercell(cell2,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
#view(cell1)
#view(cell2)

os.makedirs('./mp1434',exist_ok=True)
os.chdir('./mp1434')
#shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
cell1 = make_supercell(cell1,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
view(cell1)
cell1.write('cell1.cif')
os.makedirs('./y-friction',exist_ok=True)
os.chdir('./y-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell1,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer1,vdirec='y')
#mk_nvt_input_deepmd_friction(cell1,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='y')
os.chdir('../')
os.makedirs('./x-friction',exist_ok=True)
os.chdir('./x-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell1,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer1,vdirec='x')
#mk_nvt_input_deepmd_friction(cell1,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='x')
os.chdir('../')
os.chdir('../')

cell2 = read('./MoS2_mp2815_slab_1_ortho.cfg')

os.makedirs('./mp2815',exist_ok=True)
os.chdir('./mp2815')
cell2 = make_supercell(cell2,([nx,0,0],[0,ny,0],[0,0,1]),wrap=True)
view(cell2)
cell2.write('cell2.cif')
os.makedirs('./y-friction',exist_ok=True)
os.chdir('./y-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell2,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer2,vdirec='y')
#mk_nvt_input_deepmd_friction(cell2,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='y')
os.chdir('../')
os.makedirs('./x-friction',exist_ok=True)
os.chdir('./x-friction')
shutil.copy(g.forcedir+'/MoS2_graph.pb','./graph.pb')
mk_nvt_input_deepmd_friction(cell2,0.0005,100,100,200000,2000000,10,300,12345,nlayer=nlayer2,vdirec='x')
#mk_nvt_input_deepmd_friction(cell2,0.0005,10,10,2000,2000,10,300,12345,nlayer=nlayer,vdirec='x')
os.chdir('../')

print(len(cell1))
print(len(cell2))

"""

"""
# ランダム吸着モデル作成
os.makedirs('./adsorp',exist_ok=True)
os.chdir('./adsorp')

#slab1_ = read('/home/A23321P/work/mySiesta/organo-metalic-complex/BOC-Mg/SiO2_slab/BLYP_Hopt/siesta.STRUCT_OUT',format='struct_out')

#complex_list = ['A','B','C']
#slab_list = [3,2,1]

mol_list = ['NaDBSA', 'DBSA', 'BYK108', 'Oleylamine', 'OleylPhosphate', 'JP502_n1', 'JP502_n2', 'JP506H_n1', 'JP506H_n2']
slab_list = ['MoS2_mp1434']

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/MoS2_mp1434.cif')
slab0 = slabgen(cell,0,0,1,10,10,1,10.0,50.0)
#os.system('atomsk slab3.cif -orthogonal-cell MoS2_mp1434_slab_3_ortho.cfg')
#view(slab[2])

datas={}
ncoord = 100
for molname in mol_list:
    os.makedirs(molname,exist_ok=True)
    os.chdir(molname)
    slab = slab0[2].copy()
    mol_tmp = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/nanoMoS_mols/'+str(molname)+'.log',format='gaussian-out')
    mol_ = mol_tmp.copy()
    lat = slab.get_cell()
    zslab=[]
    for atom in slab:
        zslab.append(float(atom.position[2]))
    zmax_slab=max(zslab)
    slab_center = [(lat[0][0]+lat[1][0])/2.0,(lat[0][1]+lat[1][1])/2.0,zmax_slab]
    mol_com = mol_.get_center_of_mass()
    mol_.translate([-mol_com[0],-mol_com[1],-mol_com[2]])
    for n in range(ncoord):
        os.makedirs('nset'+str(n+1),exist_ok=True)
        os.chdir('nset'+str(n+1))
        xtheta = random.uniform(-180,180)
        ytheta = random.uniform(-180,180)
        ztheta = random.uniform(-180,180)
        height = random.uniform(3.0,3.1)
        mol0 = mol_.copy()
        mol0.rotate(xtheta,'x')
        mol0.rotate(ytheta,'y')
        mol0.rotate(ztheta,'z')
        zmol=[]
        for atom in mol0:
            zmol.append(atom.position[2])
        zmin_mol = min(zmol)
        mol0.translate([slab_center[0],slab_center[1],(zmax_slab-zmin_mol)+height])
        mol0.set_cell(lat)
        adsorp = slab+mol0
        adsorp = sortmol(adsorp,sort_atom=False)
        #view(adsorp)
        #print(len(adsorp))
        mk_mimize_input_uff_adsorp_onMoS2(adsorp,200000,100000,100000)
        adsorp.write(molname+'_adsorp_'+str(n)+'.cif')
        os.system('mpirun -np 20 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
        f = open('log_lammps','r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if 'Potential Energy' in line:
                line = line.split()
                print(line[2])
                datas[('mol_'+molname,'nset'+str(n+1),xtheta,ytheta,ztheta)] = float(line[2])
                break
        #sys.exit()    
        os.chdir('../')
    os.chdir('../')    
df = pd.Series(datas)
df.index.names = ['mol','nset','xtheta','ytheta','ztheta']
df.to_csv('potential_datas.csv')


# ランダム吸着モデル作成 終了 
"""

# 吸着モデル構造最適化(Siesta)
os.makedirs('./adsorp',exist_ok=True)
os.chdir('./adsorp')

os.makedirs('./opt',exist_ok=True)
os.chdir('./opt')

os.makedirs('./MoS2_slab',exist_ok=True)
os.chdir('./MoS2_slab')

cell = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/MoS2_mp1434.cif')
slab0 = slabgen(cell,0,0,1,10,10,1,10.0,50.0)
slab = slab0[2].copy()
#view(slab)
#sys.exit()

com = slab.get_center_of_mass()
constraint=[]
for atom in slab:
    if atom.position[2] < 14.0:
        constraint.append(('atom',atom.index+1))
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
options['Geometry.Constraints'] = constraint
options['DFTD3'] = True
options['DFTD3.UseXCDefaults'] = True
#mk_siesta_input_scf(slab,'PBE','DZP',200.0,[1,1,1],g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
mk_siesta_input_optimize(slab,'PBE','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
slab.write('test.cif') 
#view(slab)
#sys.exit()
os.chdir('../')

adsorps={}
#adsorps['NaDBSA'] = [37,56,97,98,84,35,17,13,59,31,68,65,11,7,3,54]
#adsorps['DBSA'] = [28,72,67,53,78,59,76,98,55,60,20,16,73,21,24,51,41,70]
#adsorps['BYK108'] = [87,83,88,22,50,10,56,97,25,69,2,96,31,6,16,91,59,40,47,39,89,63,15,17,84,55,60,85,67,98,61]
#adsorps['Oleylamine'] = [67,13,89,81,46,86,77,49,6,84,54,15,100,37,16,26,52]
#adsorps['OleylPhosphate'] = [34,57,97,9,27,11,10,94,81,55,95,46,16,30,70]
#adsorps['JP502_n1'] = [88,77,28,71,89,29,20,97,42,49,96]
#adsorps['JP502_n2'] = [33,61,76,94,50,74,35,46,95,10,91,88]
#adsorps['JP506H_n1'] = [28,35,73,77,12,95,92]
#adsorps['JP506H_n2'] = [80,9,14,76,23,75,36,17,5,70]

adsorps['NaDBSA'] = [98,54]
adsorps['DBSA'] = [72,70]
adsorps['BYK108'] = [83]
adsorps['Oleylamine'] = [67,13,89,81,46,86,77,49,6,84,54,15,100,37,16,26,52]
adsorps['OleylPhosphate'] = [34,57,97,9,27,11,10,94,81,55,95,46,16,30,70]
adsorps['JP502_n1'] = [88,77,28,71,89,29,20,97,42,49,96]
adsorps['JP502_n2'] = [33,61,76,94,50,74,35,46,95,10,91,88]
adsorps['JP506H_n1'] = [28,35,73,77,12,95,92]
adsorps['JP506H_n2'] = [80,9,14,76,23,75,36,17,5,70]

for molname in adsorps:
    os.makedirs(molname,exist_ok=True)
    os.chdir(molname)
    os.makedirs('mol_opt',exist_ok=True)
    os.chdir('mol_opt')
    mol_ = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/nanoMoS_mols/'+str(molname)+'.log',format='gaussian-out')
    x=[]
    y=[]
    z=[]
    for atom in mol_:
        x.append(atom.position[0])
        y.append(atom.position[1])
        z.append(atom.position[2])
    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    zmin = min(z)
    zmax = max(z)
    x_lat = xmax-xmin+15.0
    y_lat = ymax-ymin+15.0
    z_lat = zmax-zmin+15.0
    mol_com = mol_.get_center_of_mass()
    mol_.translate([-mol_com[0],-mol_com[1],-mol_com[2]])
    mol_.translate([x_lat/2.0,y_lat/2.0,z_lat/2.0])
    lat=mol_.get_cell()
    lat[0][0] = x_lat
    lat[0][1] = 0.0
    lat[0][2] = 0.0
    lat[1][0] = 0.0
    lat[1][1] = y_lat
    lat[1][2] = 0.0
    lat[2][0] = 0.0
    lat[2][1] = 0.0
    lat[2][2] = z_lat
    mol_.set_cell(cell=lat)
    options={}
    options['DFTD3'] = True
    options['DFTD3.UseXCDefaults'] = True
    #mk_siesta_input_scf(mol_,'PBE','DZP',200.0,[1,1,1],g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
    mk_siesta_input_optimize(mol_,'PBE','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
    mol_.write('test.cif') 
    #view(mol_)
    os.chdir('../')
    os.makedirs('adsorp_opt',exist_ok=True)
    os.chdir('adsorp_opt')
    for iset in adsorps[molname]:
        os.makedirs('nset_'+str(iset),exist_ok=True)
        os.chdir('nset_'+str(iset))
        adsorp = read('/home/A23321P/work/myPython/AtomicVirtuaLab/MoS2_friction/adsorp/'+molname+'/'+molname+'_nset'+str(iset)+'.cif')
        constraint=[]
        for atom in adsorp:
            if atom.position[2] < 14.0:
                constraint.append(('atom',atom.index+1))
        options={}
        options['Slab.DipoleCorrection'] = True
        options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
        options['Geometry.Constraints'] = constraint
        options['DFTD3'] = True
        options['DFTD3.UseXCDefaults'] = True
        #mk_siesta_input_scf(adsorp,'PBE','DZP',200.0,[1,1,1],g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
        mk_siesta_input_optimize(adsorp,'PBE','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
        adsorp.write('test.cif')
        #view(adsorp)
        #sys.exit()
        os.chdir('../')
    os.chdir('../')
    os.chdir('../')         

# 吸着モデル構造最適化(Siesta) 終了
