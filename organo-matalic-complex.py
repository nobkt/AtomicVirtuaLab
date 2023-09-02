from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.siesta import mk_siesta_input_optimize
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.build import sortmol
from ase.io import read, write
from ase.visualize import view
from ase import Atoms
import os
import math

g.siesta_pot = '/home/A23321P/work/myPython/AtomicVirtuaLab/siesta_pseudo'

os.makedirs('organo-metalic-complex',exist_ok=True)
os.chdir('organo-metalic-complex')

boc_smiles = 'CC(C)(C)OC(=O)NC(CC(=O)OC1CCCCC1)C([O-])=O'

smiles2xyz(boc_smiles,'boc',True,smarts=False,userandom=True)

boc1 = read('./boc.xyz')

mg = Atoms('Mg', [[0.0, 0.0, 0.0]])

lat = [50.0,50.0,50.0]

boc1.set_cell(lat)
mg.set_cell(lat)

## aligne
id0 = 8
id1 = 19
shift = boc1[id0].position
boc1.translate([-shift[0],-shift[1],-shift[2]])

#view(boc1)

vec = boc1[id1].position - boc1[id0].position
vec_xy = [vec[0],vec[1],0]
norm_vec_xy = math.sqrt(vec_xy[0]**2+vec_xy[1]**2+vec_xy[2]**2)
theta_z = math.acos(vec_xy[0]/norm_vec_xy)*180.0/math.pi
boc1.rotate(theta_z,'-z')

#view(boc1)

vec = boc1[id1].position - boc1[id0].position
vec_xz = [vec[0],0.0,vec[2]]
norm_vec_xz = math.sqrt(vec_xz[0]**2+vec_xz[1]**2+vec_xz[2]**2)
theta_y = math.acos(vec_xz[0]/norm_vec_xz)*180.0/math.pi
boc1.rotate(theta_y,'y')

#view(boc1)

id2 = 21

vec = boc1[id2].position - boc1[id0].position
vec_yz = [0,vec[1],vec[2]]
norm_vec_yz = math.sqrt(vec_yz[0]**2+vec_yz[1]**2+vec_yz[2]**2)
theta_x = math.acos(vec_yz[1]/norm_vec_yz)*180.0/math.pi
boc1.rotate(theta_x,'-x')

#view(boc1)

shift = boc1[id1].position
boc1.translate([lat[0]/2.0-4.5,lat[1]/2.0,lat[2]/2.0])

#view(boc1)

boc2 = boc1.copy()
boc2.rotate(180.0,'y',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))
boc2.rotate(180.0,'x',center=(lat[0]/2.0,lat[1]/2.0,lat[2]/2.0))

#view(boc2)

mg.translate([lat[0]/2.0,lat[1]/2.0,lat[2]/2.0])
#boc1.translate([0.0,0.0,lat[2]/2.0])
#boc2.translate([lat[0],0.0,lat[2]/2.0])

complex = mg+boc1+boc2
complex = sortmol(complex)
view(complex)

for xc in ['BLYP','PBE','BH']:
    os.makedirs(xc+'_opt',exist_ok=True)
    os.chdir(xc+'_opt')
    mk_siesta_input_optimize(complex,xc,'DZP',300.0,[1,1,1],2000,g.siesta_pot,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized')
    os.chdir('../')
