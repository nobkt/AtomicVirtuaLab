from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol, slabgen
from AtomicVirtuaLab.lammps import mk_mimize_input_uff_adsorp,mk_opt_input_deepmd
from ase.io import read, write
from ase.visualize import view
from ase import Atom, Atoms
import os
import math
import sys
import random
import pandas as pd
import shutil

g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.forcedir = "/home/A23321P/work/myPython/AtomicVirtuaLab/lmp_potentials"


os.makedirs('organo-metalic-complex',exist_ok=True)
os.chdir('organo-metalic-complex')

"""
# 分散剤の構造最適化

lat = [50.0,50.0,50.0]

os.makedirs('mols',exist_ok=True)
os.chdir('mols')
os.makedirs('dispersant',exist_ok=True)
os.chdir('dispersant')

os.makedirs('AD3284M',exist_ok=True)
os.chdir('AD3284M')

smiles = 'CC(COCCO)OCN(COC(C)COCCO)COC(C)COCCO'

smiles2xyz(smiles,'mol_AD3284M',True,smarts=False,userandom=True)

mol_AD3284M = read('./mol_AD3284M.xyz')
mol_AD3284M.set_cell(lat)
com = mol_AD3284M.get_center_of_mass()
shift = [lat[0]/2.0-com[0],lat[1]/2.0-com[1],lat[2]/2.0-com[2]]
mol_AD3284M.translate(shift)
mk_siesta_input_optimize(mol_AD3284M,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
view(mol_AD3284M)

os.chdir('../')

os.makedirs('BYK191',exist_ok=True)
os.chdir('BYK191')

smiles = '[H]CC(CC(CC(CC(C(\C=C/CCCCCC)C(CCCCCCC(O)=O)C(CCCCCCCCC)C([H])CCCCCCC(O)=O)C(=O)OCCCC)C(=O)OCCN(C)C)C(O)=O)C(=O)OCCOCCOC'

smiles2xyz(smiles,'mol_BYK191',True,smarts=False,userandom=True)

mol_BYK191 = read('./mol_BYK191.xyz')
mol_BYK191.set_cell(lat)
com = mol_BYK191.get_center_of_mass()
shift = [lat[0]/2.0-com[0],lat[1]/2.0-com[1],lat[2]/2.0-com[2]]
mol_BYK191.translate(shift)
mk_siesta_input_optimize(mol_BYK191,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
view(mol_BYK191)

os.chdir('../')

os.makedirs('S71000',exist_ok=True)
os.chdir('S71000')

smiles = 'COCCOC(C)COC(=O)CCC(=O)NCCNCCN(CCNCCNC(=O)CCC(=O)OCC(C)OCCOC)C(=O)CCC(=O)OCC(C)OCCOC'

smiles2xyz(smiles,'mol_S71000',True,smarts=False,userandom=True)

mol_S71000 = read('./mol_S71000.xyz')
mol_S71000.set_cell(lat)
com = mol_S71000.get_center_of_mass()
shift = [lat[0]/2.0-com[0],lat[1]/2.0-com[1],lat[2]/2.0-com[2]]
mol_S71000.translate(shift)
mk_siesta_input_optimize(mol_S71000,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
view(mol_S71000)

os.chdir('../')

os.makedirs('S75000',exist_ok=True)
os.chdir('S75000')

smiles = 'CCCCCCCCCCC(=O)OCCCCCCC(=O)NCCNCCN(CCNCCNC(=O)CCCCCCOC(=O)CCCCCCCCCC)C(=O)CCCCCCOC(=O)CCCCCCCCCC'

smiles2xyz(smiles,'mol_S75000',True,smarts=False,userandom=True)

mol_S75000 = read('./mol_S75000.xyz')
mol_S75000.set_cell(lat)
com = mol_S75000.get_center_of_mass()
shift = [lat[0]/2.0-com[0],lat[1]/2.0-com[1],lat[2]/2.0-com[2]]
mol_S75000.translate(shift)
mk_siesta_input_optimize(mol_S75000,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
view(mol_S75000)

os.chdir('../')

# 分散剤の構造最適化 終了
"""

"""
# optimize molA-Mg
os.makedirs('mols',exist_ok=True)
os.chdir('mols')
os.makedirs('molA',exist_ok=True)
os.chdir('molA')

molA_smiles = 'CCCCCC(CC)C([O-])=O'

smiles2xyz(molA_smiles,'molA',True,smarts=False,userandom=True)

molA1 = read('./molA.xyz')
#view(molA1)
#sys.exit()

mg = Atoms('Mg', [[0.0, 0.0, 0.0]])

lat = [50.0,50.0,50.0]

molA1.set_cell(lat)
mg.set_cell(lat)
#view(molA1)

## aligne
id0 = 8
id1 = 5
shift = molA1[id0].position
molA1.translate([-shift[0],-shift[1],-shift[2]])

#view(molA1)
#sys.exit()

vec = molA1[id1].position -molA1[id0].position
vec_xy = [vec[0],vec[1],0]
norm_vec_xy = math.sqrt(vec_xy[0]**2+vec_xy[1]**2+vec_xy[2]**2)
theta_z = math.acos(vec_xy[0]/norm_vec_xy)*180.0/math.pi
molA1.rotate(theta_z,'-z')

#view(molA1)
#sys.exit()

vec = molA1[id1].position - molA1[id0].position
vec_xz = [vec[0],0.0,vec[2]]
norm_vec_xz = math.sqrt(vec_xz[0]**2+vec_xz[1]**2+vec_xz[2]**2)
theta_y = math.acos(vec_xz[0]/norm_vec_xz)*180.0/math.pi
molA1.rotate(theta_y,'y')

#view(molA1)
#sys.exit()

id2 = 9

vec = molA1[id2].position - molA1[id0].position
vec_yz = [0,vec[1],vec[2]]
norm_vec_yz = math.sqrt(vec_yz[0]**2+vec_yz[1]**2+vec_yz[2]**2)
theta_x = math.acos(vec_yz[1]/norm_vec_yz)*180.0/math.pi
molA1.rotate(theta_x,'-x')

#view(molA1)
#sys.exit()

shift = molA1[id1].position
molA1.translate([lat[0]/2.0+4.5,lat[1]/2.0,lat[2]/2.0])

#view(molA1)
#sys.exit()

molA2 = molA1.copy()
molA2.rotate(180.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molA2.rotate(0.0,'x',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))

#view(molA2)

mg.translate([lat[0]/2.0,lat[1]/2.0,lat[2]/2.0])

complex = mg+molA1+molA2
#complex = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/molA_Quench.mol')
complex.set_cell(lat)
complex = sortmol(complex)
#view(complex)
#sys.exit()
#shift = [lat[0]/2.0-complex[52].position[0],lat[1]/2.0-complex[52].position[1],lat[2]/2.0-complex[52].position[2]]
#complex.translate(shift)
view(complex)
#complex = sortmol(complex)
#view(complex)
#sys.exit()

for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    os.makedirs('molA_opt',exist_ok=True)
    os.chdir('molA_opt')
    complex.write('test.xyz',format='xyz')
    mk_siesta_input_optimize(complex,xc,'DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')
    os.chdir('../')

os.chdir('../')
os.chdir('../')
# optimize molA-Mg 終了


# optimize molB-Mg
os.makedirs('mols',exist_ok=True)
os.chdir('mols')
os.makedirs('molB',exist_ok=True)
os.chdir('molB')

molB_smiles = 'CC(C)(C)OC(=O)NC(CC(=O)OC1CCCCC1)C([O-])=O'

smiles2xyz(molB_smiles,'molB',True,smarts=False,userandom=True)

molB1 = read('./molB.xyz')

mg = Atoms('Mg', [[0.0, 0.0, 0.0]])

lat = [50.0,50.0,50.0]

molB1.set_cell(lat)
mg.set_cell(lat)

## aligne
id0 = 8
id1 = 19
shift = molB1[id0].position
molB1.translate([-shift[0],-shift[1],-shift[2]])

#view(molB1)

vec = molB1[id1].position -molB1[id0].position
vec_xy = [vec[0],vec[1],0]
norm_vec_xy = math.sqrt(vec_xy[0]**2+vec_xy[1]**2+vec_xy[2]**2)
theta_z = math.acos(vec_xy[0]/norm_vec_xy)*180.0/math.pi
molB1.rotate(theta_z,'-z')

#view(molB1)
#sys.exit()

vec = molB1[id1].position - molB1[id0].position
vec_xz = [vec[0],0.0,vec[2]]
norm_vec_xz = math.sqrt(vec_xz[0]**2+vec_xz[1]**2+vec_xz[2]**2)
theta_y = math.acos(vec_xz[0]/norm_vec_xz)*180.0/math.pi
molB1.rotate(theta_y,'-y')

#view(molB1)
#sys.exit()

id2 = 21

vec = molB1[id2].position - molB1[id0].position
vec_yz = [0,vec[1],vec[2]]
norm_vec_yz = math.sqrt(vec_yz[0]**2+vec_yz[1]**2+vec_yz[2]**2)
theta_x = math.acos(vec_yz[1]/norm_vec_yz)*180.0/math.pi
molB1.rotate(theta_x,'-x')

#view(molB1)
#sys.exit()

shift = molB1[id1].position
molB1.translate([lat[0]/2.0-4.5,lat[1]/2.0,lat[2]/2.0])

#view(molB1)

molB2 = molB1.copy()
molB2.rotate(180.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molB2.rotate(180.0,'x',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))

#view(molB2)

mg.translate([lat[0]/2.0,lat[1]/2.0,lat[2]/2.0])

complex = mg+molB1+molB2
#complex = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/molB_Quench.mol')
complex.set_cell(lat)
complex = sortmol(complex)
#view(complex)
#sys.exit()
shift = [lat[0]/2.0-complex[78].position[0],lat[1]/2.0-complex[78].position[1],lat[2]/2.0-complex[78].position[2]]
complex.translate(shift)
view(complex)
#complex = sortmol(complex)
#view(complex)
complex.write('test.xyz',format='xyz')
sys.exit()

for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    os.makedirs('MolB_opt',exist_ok=True)
    os.chdir('MolB_opt')
    complex.write('test.xyz',format='xyz')
    mk_siesta_input_optimize(complex,xc,'DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')
    os.chdir('../')

os.chdir('../')
os.chdir('../')
# optimize MolB-Mg 終了


# optimize molC-Mg
os.makedirs('mols',exist_ok=True)
os.chdir('mols')
os.makedirs('molC',exist_ok=True)
os.chdir('molC')

molC_smiles = 'CCCCCC(CC)C([O-])=O'
molC0_smiles = 'C1=CN=C(C=C1)C1=NC=CC=C1'

smiles2xyz(molC_smiles,'molC',True,smarts=False,userandom=True)
smiles2xyz(molC0_smiles,'molC0',True,smarts=False,userandom=True)

molC1 = read('./molC.xyz')
molC0 = read('./molC0.xyz')
#view(molC1)
#sys.exit()

mg = Atoms('Mg', [[0.0, 0.0, 0.0]])

lat = [50.0,50.0,50.0]

molC1.set_cell(lat)
molC0.set_cell(lat)
mg.set_cell(lat)
#view(molC1)
#view(molC0)
#sys.exit()

## aligne
id0 = 8
id1 = 5
shift = molC1[id0].position
molC1.translate([-shift[0],-shift[1],-shift[2]])

jd0 = 2
jd1 = 7
shift = molC0[jd0].position
molC0.translate([-shift[0],-shift[1],-shift[2]])

#view(molC0)
#sys.exit()

vec = molC1[id1].position -molC1[id0].position
vec_xy = [vec[0],vec[1],0]
norm_vec_xy = math.sqrt(vec_xy[0]**2+vec_xy[1]**2+vec_xy[2]**2)
theta_z = math.acos(vec_xy[0]/norm_vec_xy)*180.0/math.pi
molC1.rotate(theta_z,'-z')

vec = molC0[jd1].position -molC0[jd0].position
vec_xy = [vec[0],vec[1],0]
norm_vec_xy = math.sqrt(vec_xy[0]**2+vec_xy[1]**2+vec_xy[2]**2)
theta_z = math.acos(vec_xy[0]/norm_vec_xy)*180.0/math.pi
molC0.rotate(theta_z,'z')
molC0.rotate(270.0,'z')
molC0.rotate(90.0,'-y')

#view(molC0)
#sys.exit()

vec = molC1[id1].position - molC1[id0].position
vec_xz = [vec[0],0.0,vec[2]]
norm_vec_xz = math.sqrt(vec_xz[0]**2+vec_xz[1]**2+vec_xz[2]**2)
theta_y = math.acos(vec_xz[0]/norm_vec_xz)*180.0/math.pi
molC1.rotate(theta_y,'y')

#view(molC1)
#sys.exit()

id2 = 9

vec = molC1[id2].position - molC1[id0].position
vec_yz = [0,vec[1],vec[2]]
norm_vec_yz = math.sqrt(vec_yz[0]**2+vec_yz[1]**2+vec_yz[2]**2)
theta_x = math.acos(vec_yz[1]/norm_vec_yz)*180.0/math.pi
molC1.rotate(theta_x,'-x')

#view(molC1)
#sys.exit()

shift = molC1[id1].position
molC1.translate([lat[0]/2.0+4.5,lat[1]/2.0,lat[2]/2.0])

shift = molC0[jd1].position
molC0.translate([lat[0]/2.0+1.5,lat[1]/2.0,lat[2]/2.0])

#view(molC0)
#sys.exit()

molC2 = molC1.copy()
molC2.rotate(0.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molC2.rotate(90.0,'x',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molC1.rotate(75.0,'z',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molC1.rotate(45.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molC0.rotate(-60.0,'z',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
molC0.rotate(-45.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))


#view(molC2)

mg.translate([lat[0]/2.0,lat[1]/2.0,lat[2]/2.0])

complex = mg+molC1+molC2+molC0
#view(complex)
complex = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/molC_Quench.mol')
complex.set_cell(lat)
complex = sortmol(complex)
shift = [lat[0]/2.0-complex[70].position[0],lat[1]/2.0-complex[70].position[1],lat[2]/2.0-complex[70].position[2]]
complex.translate(shift)
view(complex)


for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    os.makedirs('molC_opt',exist_ok=True)
    os.chdir('molC_opt')
    complex.write('test.xyz',format='xyz')
    mk_siesta_input_optimize(complex,xc,'DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')
    os.chdir('../')

os.chdir('../')
os.chdir('../')
# optimize molC-Mg 終了
#sys.exit()
"""

"""
# SiO2 001 slab
os.makedirs('./SiO2_slab',exist_ok=True)
os.chdir('./SiO2_slab')

cell = read(g.cifdir+'/SiO2_quartz_primitive.cif',primitive_cell=True,subtrans_included=False)
#view(cell)

slab1 = slabgen(cell,1,0,-1,5,3,2,5.0,25.0)
slab2 = slabgen(cell,1,0,-1,6,4,2,5.0,25.0)
slab3 = slabgen(cell,1,0,-1,7,5,2,5.0,25.0)
slab4 = slabgen(cell,1,0,-1,8,6,2,5.0,25.0)
#view(slab1[0])
#view(slab2[0])
#view(slab3[0])
#view(slab4[0])

slab1_ = slab1[0].copy()
slab2_ = slab2[0].copy()
slab3_ = slab3[0].copy()
slab4_ = slab4[0].copy()

rOH = 0.9
z0 = 10.0
for atom in slab1[0]:
    if atom.symbol == 'O' and atom.position[2] > z0:
        lowH = Atom('H',(atom.position[0],atom.position[1],atom.position[2]+rOH))
        slab1_.append(lowH)

for atom in slab2[0]:
    if atom.symbol == 'O' and atom.position[2] > z0:
        lowH = Atom('H',(atom.position[0],atom.position[1],atom.position[2]+rOH))
        slab2_.append(lowH)

for atom in slab3[0]:
    if atom.symbol == 'O' and atom.position[2] > z0:
        lowH = Atom('H',(atom.position[0],atom.position[1],atom.position[2]+rOH))
        slab3_.append(lowH)

for atom in slab4[0]:
    if atom.symbol == 'O' and atom.position[2] > z0:
        lowH = Atom('H',(atom.position[0],atom.position[1],atom.position[2]+rOH))
        slab4_.append(lowH)
    #slab4_ = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/SiO2_101_model4.cif')

#view(slab1_)
#view(slab2_)
#view(slab3_)
#view(slab4_)

for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    i = 1
    for slabs in [slab1_,slab2_,slab3_,slab4_]:
        os.makedirs('slab'+str(i)+'_opt',exist_ok=True)
        os.chdir('slab'+str(i)+'_opt')
        com = slabs.get_center_of_mass()
        constraint=[]
        for atom in slabs:
            if atom.position[2] < 5.3:
                constraint.append(('atom',atom.index+1))
        options={}
        options['Slab.DipoleCorrection'] = True
        options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
        options['Geometry.Constraints'] = constraint
        mk_siesta_input_optimize(slabs,xc,'DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
        slabs.write('test.cif')
        os.chdir('../')
        i=i+1
    os.chdir('../')
os.chdir('../')
#sys.exit()
# SiO2 001 slab 終了
sys.exit()
"""

"""
# ランダム吸着モデル作成
os.makedirs('./adsorp1',exist_ok=True)
os.chdir('./adsorp1')

#slab1_ = read('/home/A23321P/work/mySiesta/organo-metalic-complex/BOC-Mg/SiO2_slab/BLYP_Hopt/siesta.STRUCT_OUT',format='struct_out')

#complex_list = ['A','B','C']
#slab_list = [3,2,1]

complex_list = ['AD3284M','BYK191','S71000','S75000']
slab_list = [3]

datas={}
ncoord = 100
for nslab in slab_list:
    for ncomplex in complex_list:
        os.makedirs('slab_'+str(nslab)+'-mol'+ncomplex,exist_ok=True)
        os.chdir('slab_'+str(nslab)+'-mol'+ncomplex)
        slab = read('/home/A23321P/work/mySiesta/organo-metalic-complex/SiO2_slab/BLYP/slab'+str(nslab)+'_opt/siesta.STRUCT_OUT',format='struct_out')
        #mol = read('/home/A23321P/work/mySiesta/organo-metalic-complex/mols/mol'+ncomplex+'/BLYP/mol'+ncomplex+'_opt/siesta.xyz',format='xyz')
        mol = read('/home/A23321P/work/mySiesta/organo-metalic-complex/mols/dispersant/'+ncomplex+'/siesta.xyz',format='xyz')
        lat = slab.get_cell()
        zslab=[]
        for atom in slab:
            zslab.append(float(atom.position[2]))
        zmax_slab=max(zslab)
        slab_center = [(lat[0][0]+lat[1][0])/2.0,(lat[0][1]+lat[1][1])/2.0,zmax_slab]
        mol_com = mol.get_center_of_mass()
        mol.translate([-mol_com[0],-mol_com[1],-mol_com[2]])
        for n in range(ncoord):
            os.makedirs('nset'+str(n+1),exist_ok=True)
            os.chdir('nset'+str(n+1))
            xtheta = random.uniform(-180,180)
            ytheta = random.uniform(-180,180)
            ztheta = random.uniform(-180,180)
            height = random.uniform(2.0,2.1)
            mol0 = mol.copy()
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
            mk_mimize_input_uff_adsorp(adsorp,200000,100000,100000)
            os.system('mpirun -np 20 lmp -in lammps.lmp 1> log_lammps 2> err_lammps')
            f = open('log_lammps','r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Potential Energy' in line:
                    line = line.split()
                    print(line[2])
                    datas[('slab'+str(nslab),'mol_'+ncomplex,'nset'+str(n+1),xtheta,ytheta,ztheta)] = float(line[2])
                    break
            #view(slab+mol0)
            #sys.exit()
            os.chdir('../')
        os.chdir('../')
df = pd.Series(datas)
df.index.names = ['slab','mol','nset','xtheta','ytheta','ztheta']
df.to_csv('potential_datas.csv')

# ランダム吸着モデル作成 終了        
"""

"""        
# ランダム吸着モデル構造最適化
os.makedirs('./adsorp1',exist_ok=True)
os.chdir('./adsorp1')

os.makedirs('./opt',exist_ok=True)
os.chdir('./opt')

complex_list = ['A','B','C']
slab_list = [1,2,3]    
ncoord = 100

#nslab0=1
#ncomp0='A'
#opt_list=[61,34,92,6,81,26,47,35,83,10]
#Z_of_type={1:6, 2:1, 3:12, 4:8, 5:14}

#ncomp0='B'
#opt_list=[99,50,46,63,80,37,36,69,67,33]
#Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

#ncomp0='C'
#opt_list=[93, 36, 35, 18, 25, 13, 1, 76, 30, 23]
#Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

#nslab0=2
#ncomp0='A'
#opt_list=[89,71,66,97,39,99,19,93,29,92]
#Z_of_type={1:6, 2:1, 3:12, 4:8, 5:14}

#ncomp0='B'
#opt_list=[70,49,85,14,74,55,32,81,25,20]
#Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

#ncomp0='C'
#opt_list=[51,99,33,13,91,27,64,18,63,32]
#Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

nslab0=3
#ncomp0='A'
#opt_list=[49,79,77,63,9,2,91,16,81,85,42,6,35,82,3,5,54,10,43,15,74]
#Z_of_type={1:6, 2:1, 3:12, 4:8, 5:14}

#ncomp0='B'
#opt_list=[89,59,33,81,5,43,60,77,95,6,65,3,51,40,1,99,93,25,91,79,28,38,82,13,67,34,71,54,68,19,55,41,22,14,49,78,17,73,42]
#Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

ncomp0='C'
opt_list=[56,31,34,92,86,4,30,51,41,91,1,52,65,63,9,10,60,16,55,58,49,97,5,64,19,90,81,96,87,61,24,77,36,54,80,98,21,22,50]
Z_of_type={1:6, 2:1, 3:12, 4:7, 5:8, 6:14}

for nslab in slab_list:
    for ncomp in complex_list:
        for n in range(ncoord):
            if nslab == nslab0 and ncomp == ncomp0 and (n+1) in opt_list:
                print(n+1)
                os.makedirs('slab_'+str(nslab)+'-mol'+ncomp,exist_ok=True)
                os.chdir('slab_'+str(nslab)+'-mol'+ncomp)
                os.makedirs('nset'+str(n+1),exist_ok=True)
                os.chdir('nset'+str(n+1))
                dir_ = '/home/A23321P/work/myPython/AtomicVirtuaLab/organo-metalic-complex/adsorp1/slab_'+str(nslab)+'-mol'+ncomp+'/nset'+str(n+1)
                adsorp = read(dir_+'/result.data',format='lammps-data',Z_of_type=Z_of_type)
                #view(adsorp)
                ii = 0
                sio2 = adsorp.copy()
                del sio2[[atom.index for atom in sio2 if sio2.arrays["mol-id"][atom.index] != 1]]
                com = sio2.get_center_of_mass()
                constraint=[]
                for atom in adsorp:
                    if atom.position[2] < 5.3:
                        constraint.append(('atom',atom.index+1))
                options={}
                options['Slab.DipoleCorrection'] = True
                options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
                options['Geometry.Constraints'] = constraint
                mk_siesta_input_optimize(adsorp,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
                adsorp.write('test.cif')                
                #sys.exit()
                os.chdir('../')
                os.chdir('../')

# ランダム吸着モデル構造最適化 終了
"""

"""
# OC20 ランダム吸着モデル構造最適化
os.makedirs('./adsorp1',exist_ok=True)
os.chdir('./adsorp1')

os.makedirs('./oc_opt',exist_ok=True)
os.chdir('./oc_opt')

complex_list = ['AD3284M','BYK191','S71000','S75000']
slab_list = [3]    
ncoord = 100

Z_of_type={1:6,2:1,3:7,4:8,5:14}

for nslab in slab_list:
    for ncomp in complex_list:
        for n in range(ncoord):
            os.makedirs('slab_'+str(nslab)+'-'+ncomp,exist_ok=True)
            os.chdir('slab_'+str(nslab)+'-'+ncomp)
            os.makedirs('nset'+str(n+1),exist_ok=True)
            os.chdir('nset'+str(n+1))
            dir_ = '/home/A23321P/work/myPython/AtomicVirtuaLab/organo-metalic-complex/adsorp1/slab_'+str(nslab)+'-mol'+ncomp+'/nset'+str(n+1)
            shutil.copy(g.forcedir+'/OC_10M.pb','./graph.pb')
            adsorp = read(dir_+'/result.data',format='lammps-data',Z_of_type=Z_of_type)
            #view(adsorp)
            mk_opt_input_deepmd(adsorp,1,1,mol=True,fixlay=5.3)
            adsorp.write('test.cif')                
            #sys.exit()
            os.chdir('../')
            os.chdir('../')

# OC20 ランダム吸着モデル構造最適化 終了
"""

"""
# ランダム吸着モデル構造最適化
os.makedirs('./adsorp1',exist_ok=True)
os.chdir('./adsorp1')

os.makedirs('./siesta_opt',exist_ok=True)
os.chdir('./siesta_opt')

complex_list = ['AD3284M','BYK191','S71000','S75000']
slab_list = [3]    
ncoord = 100

Z_of_type={1:6,2:1,3:7,4:8,5:14}

opt_list={}
opt_list['AD3284M'] = [48]
opt_list['BYK191'] = [4]
opt_list['S71000'] = [25]
opt_list['S75000'] = [55]

for nslab in slab_list:
    for ncomp in complex_list:
        for n in range(ncoord):
            if n+1 not in opt_list[ncomp]:
                continue
            os.makedirs('slab_'+str(nslab)+'-'+ncomp,exist_ok=True)
            os.chdir('slab_'+str(nslab)+'-'+ncomp)
            os.makedirs('nset'+str(n+1),exist_ok=True)
            os.chdir('nset'+str(n+1))
            dir_ = '/home/A23321P/work/myPython/AtomicVirtuaLab/organo-metalic-complex/adsorp1_init/slab_'+str(nslab)+'-mol'+ncomp+'/nset'+str(n+1)
            adsorp = read(dir_+'/result.data',format='lammps-data',Z_of_type=Z_of_type)
            view(adsorp)
            sio2 = adsorp.copy()
            del sio2[[atom.index for atom in sio2 if sio2.arrays["mol-id"][atom.index] != 1]]
            com = sio2.get_center_of_mass()
            constraint=[]
            for atom in adsorp:
                if atom.position[2] < 5.3:
                    constraint.append(('atom',atom.index+1))
            options={}
            options['Slab.DipoleCorrection'] = True
            options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
            options['Geometry.Constraints'] = constraint
            mk_siesta_input_optimize(adsorp,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
            adsorp.write('test.cif')             
            #sys.exit()
            os.chdir('../')
            os.chdir('../')

# ランダム吸着モデル構造最適化 終了
"""

# テスト
os.makedirs('test',exist_ok=True)
os.chdir('test')

adsorp = read(g.cifdir+'/test/test.cif')
sio2 = adsorp.copy()
del sio2[[atom.index for atom in sio2 if atom.position[2] >= 11.571]]
com = sio2.get_center_of_mass()
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
mk_siesta_input_optimize(adsorp,'BLYP','DZP',200.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
view(adsorp)