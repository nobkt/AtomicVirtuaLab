from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol, slabgen
from ase.io import read, write
from ase.visualize import view
from ase import Atom, Atoms
import os
import math
import sys

g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'
g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'


os.makedirs('organo-metalic-complex',exist_ok=True)
os.chdir('organo-metalic-complex')

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
complex = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/molA_Quench.mol')
complex.set_cell(lat)
complex = sortmol(complex)
#view(complex)
#sys.exit()
shift = [lat[0]/2.0-complex[52].position[0],lat[1]/2.0-complex[52].position[1],lat[2]/2.0-complex[52].position[2]]
complex.translate(shift)
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
    mk_siesta_input_optimize(complex,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
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
complex = read('/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/molB_Quench.mol')
complex.set_cell(lat)
complex = sortmol(complex)
#view(complex)
#sys.exit()
shift = [lat[0]/2.0-complex[78].position[0],lat[1]/2.0-complex[78].position[1],lat[2]/2.0-complex[78].position[2]]
complex.translate(shift)
view(complex)
#complex = sortmol(complex)
#view(complex)
#sys.exit()

for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    os.makedirs('MolB_opt',exist_ok=True)
    os.chdir('MolB_opt')
    complex.write('test.xyz',format='xyz')
    mk_siesta_input_optimize(complex,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
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
#sys.exit()

for xc in ['BLYP']:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    os.makedirs('molC_opt',exist_ok=True)
    os.chdir('molC_opt')
    complex.write('test.xyz',format='xyz')
    mk_siesta_input_optimize(complex,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')
    os.chdir('../')

os.chdir('../')
os.chdir('../')
# optimize molC-Mg 終了
#sys.exit()

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
        mk_siesta_input_optimize(slabs,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
        slabs.write('test.cif')
        os.chdir('../')
        i=i+1
    os.chdir('../')
os.chdir('../')
#sys.exit()
# SiO2 001 slab 終了
sys.exit()
# 吸着モデル作成
os.makedirs('./adsorp',exist_ok=True)
os.chdir('./adsorp')

#slab1_ = read('/home/A23321P/work/mySiesta/organo-metalic-complex/BOC-Mg/SiO2_slab/BLYP_Hopt/siesta.STRUCT_OUT',format='struct_out')

complex1 = complex.copy()
complex2 = complex.copy()
complex3 = complex.copy()
complex4 = complex.copy()

slab1_lat = slab1_.get_cell()
complex1.set_cell(slab1_lat)

slab2_lat = slab2_.get_cell()
complex2.set_cell(slab2_lat)

slab3_lat = slab3_.get_cell()
complex3.set_cell(slab3_lat)

slab4_lat = slab4_.get_cell()
complex4.set_cell(slab4_lat)

#complex1.rotate(90.0,'y',center=(25,25,25))
complex1.rotate(180.0,'x',center=(25,25,25))
complex1.rotate(150.0,'z',center=(25,25,25))

#complex1.rotate(90.0,'y',center=(25,25,25))
complex2.rotate(180.0,'x',center=(25,25,25))
complex2.rotate(150.0,'z',center=(25,25,25))

#complex1.rotate(90.0,'y',center=(25,25,25))
complex3.rotate(180.0,'x',center=(25,25,25))
complex3.rotate(150.0,'z',center=(25,25,25))

#complex1.rotate(90.0,'y',center=(25,25,25))
complex4.rotate(180.0,'x',center=(25,25,25))
complex4.rotate(150.0,'z',center=(25,25,25))

shift1 = slab1_[140].position - complex1[78].position
shift2 = slab2_[264].position - complex2[78].position
shift3 = slab3_[320].position - complex3[78].position
shift4 = slab4_[498].position - complex4[78].position

#adsorp 1
shift1[2] = shift1[2]+4.5
shift2[2] = shift2[2]+4.5
shift3[2] = shift3[2]+4.5
shift4[2] = shift4[2]+4.5

#adsorp 2
#shift[2] = shift[2]+11.0

#adsorp 3
#shift[2] = shift[2]+15.5

complex1.translate(shift1)
complex2.translate(shift2)
complex3.translate(shift3)
complex4.translate(shift4)

adsorp1 = slab1_+complex1
os.makedirs('./adsorp1_opt',exist_ok=True)
os.chdir('./adsorp1_opt')
com = slab1_.get_center_of_mass()
constraint=[]
for atom in adsorp1:
    if atom.position[2] < 5.3:
        constraint.append(('atom',atom.index+1))
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
options['Geometry.Constraints'] = constraint
mk_siesta_input_optimize(adsorp1,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
view(adsorp1)
os.chdir('../')

adsorp2 = slab2_+complex2
os.makedirs('./adsorp2_opt',exist_ok=True)
os.chdir('./adsorp2_opt')
com = slab2_.get_center_of_mass()
constraint=[]
for atom in adsorp2:
    if atom.position[2] < 5.3:
        constraint.append(('atom',atom.index+1))
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
options['Geometry.Constraints'] = constraint
mk_siesta_input_optimize(adsorp2,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
view(adsorp2)
os.chdir('../')

adsorp3 = slab3_+complex3
os.makedirs('./adsorp3_opt',exist_ok=True)
os.chdir('./adsorp3_opt')
com = slab3_.get_center_of_mass()
constraint=[]
for atom in adsorp3:
    if atom.position[2] < 5.3:
        constraint.append(('atom',atom.index+1))
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
options['Geometry.Constraints'] = constraint
mk_siesta_input_optimize(adsorp3,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
view(adsorp3)
os.chdir('../')

adsorp4 = slab4_+complex4
os.makedirs('./adsorp4_opt',exist_ok=True)
os.chdir('./adsorp4_opt')
com = slab4_.get_center_of_mass()
constraint=[]
for atom in adsorp4:
    if atom.position[2] < 5.3:
        constraint.append(('atom',atom.index+1))
options={}
options['Slab.DipoleCorrection'] = True
options['Slab.DipoleCorrection.Origin'] = [' '+str(com[0])+' '+str(com[1])+' '+str(com[2])+' Ang']
options['Geometry.Constraints'] = constraint
mk_siesta_input_optimize(adsorp3,xc,'DZP',50.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,options=options,spin='non-polarized')
view(adsorp4)
adsorp4.write('test.cif')
os.chdir('../')

os.chdir('../')

# 吸着モデル作成 終了