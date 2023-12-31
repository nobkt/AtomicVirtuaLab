from ase.io import read
from ase.visualize import view
from ase.calculators.vasp import Vasp
import AtomicVirtuaLab.globalv as g
import os

g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

os.makedirs('UO2_PuO2_AmO2',exist_ok=True)
os.chdir('UO2_PuO2_AmO2')

"""
# UO2 cell opt
os.makedirs('UO2',exist_ok=True)
os.chdir('UO2')
cell = read(g.cifs+'/UO2.cif')
zn=14.0
#view(cell)

startz = 0.0
dz = 2.727
nz = 2

magmom=[]
for atom in cell:
    if atom.symbol == 'U':
        for n in range(nz):
            d0 = startz+dz*float(n)
            d1 = d0 + dz
            if n%2 == 0:
                direc = 1.0
            else:
                direc = -1.0
            if d0 <= atom.position[2] and atom.position[2] < d1:
                magmom.append(direc)
    else:
        magmom.append(0.0)
#print(magmom)

cell.set_initial_magnetic_moments(magmom)

lencut = [400,450,500,550,600,650]
lkpts = [4,5,6]
lplusU = [1.0,2.0,3.0,4.0,5.0]

os.makedirs('optcell',exist_ok=True)
os.chdir('optcell')

for encut in lencut:
    os.makedirs('encut_'+str(encut),exist_ok=True)
    os.chdir('encut_'+str(encut))
    for kk in lkpts:
        os.makedirs('kpt_'+str(kk),exist_ok=True)
        os.chdir('kpt_'+str(kk))
        for plusU in lplusU:
            os.makedirs('plusU_'+str(plusU),exist_ok=True)
            os.chdir('plusU_'+str(plusU))
            ldau_luj={'U': {'L': 3, 'U': plusU, 'J': 0}}
            kpts=[kk,kk,kk]
            if plusU != 0.0:
                Vasp(xc='scan',
                     ismear=0,
                     ncore=4,
                     lreal='Auto',
                     algo='Normal',
                     prec='Accurate',
                     lasph=True,
                     isym=0,
                     ibrion=2,
                     isif=3,
                     nsw=10000,
                     nelm=10000,
                     amix=0.1,
                     encut=encut,
                     ldau_luj=ldau_luj,
                     lmaxmix=6,
                     ldautype=1,
                     kpts=kpts,
                     gamma=True).write_input(cell)
            else:
                Vasp(xc='scan',
                     ismear=0,
                     ncore=4,
                     lreal='Auto',
                     algo='Normal',
                     prec='Accurate',
                     lasph=True,
                     isym=0,
                     ibrion=2,
                     isif=3,
                     nsw=10000,
                     nelm=10000,
                     amix=0.1,
                     encut=encut,
                     kpts=kpts,
                     gamma=True).write_input(cell)
            os.chdir('../')
        os.chdir('../')
    os.chdir('..')
os.chdir('../')
"""

# UO2 elastic constant
os.makedirs('UO2',exist_ok=True)
os.chdir('UO2')

os.makedirs('elastic',exist_ok=True)
os.chdir('elastic')

cell = read(g.cifs+'/UO2.cif')
#zn=14.0
#view(cell)

startz = 0.0
dz = 2.727
nz = 2

magmom=[]
for atom in cell:
    if atom.symbol == 'U':
        for n in range(nz):
            d0 = startz+dz*float(n)
            d1 = d0 + dz
            if n%2 == 0:
                direc = 1.0
            else:
                direc = -1.0
            if d0 <= atom.position[2] and atom.position[2] < d1:
                magmom.append(direc)
    else:
        magmom.append(0.0)
#print(magmom)

lencut = [400,450,500,550,600,650]
lkpts = [4,5,6]
lplusU = [1.0,2.0,3.0,4.0,5.0]

for encut in lencut:
    for kpt in lkpts:
        for plusU in lplusU:
            cell = read('/home/A23321P/work/myVASP/UO2_PuO2_AmO2/UO2/optcell/encut_'+str(encut)+'/kpt_'+str(kpt)+'/plusU_'+str(plusU)+'/CONTCAR')
            cell.set_initial_magnetic_moments(magmom)
            os.makedirs('encut_'+str(encut),exist_ok=True)
            os.chdir('encut_'+str(encut))
            os.makedirs('kpt_'+str(kpt),exist_ok=True)
            os.chdir('kpt_'+str(kpt))
            os.makedirs('plusU_'+str(plusU),exist_ok=True)
            os.chdir('plusU_'+str(plusU))
            ldau_luj={'U': {'L': 3, 'U': plusU, 'J': 0}}
            kpts=[kpt,kpt,kpt]
            if plusU != 0.0:
                Vasp(xc='scan',
                     ismear=0,
                     lreal='Auto',
                     algo='Normal',
                     prec='Accurate',
                     lasph=True,
                     isym=0,
                     ibrion=6,
                     isif=3,
                     nelm=10000,
                     amix=0.1,
                     encut=encut,
                     ldau_luj=ldau_luj,
                     lmaxmix=6,
                     ldautype=1,
                     kpts=kpts,
                     gamma=True).write_input(cell)
            else:
                Vasp(xc='scan',
                     ismear=0,
                     lreal='Auto',
                     algo='Normal',
                     prec='Accurate',
                     lasph=True,
                     isym=0,
                     ibrion=5,
                     isif=3,
                     nelm=10000,
                     amix=0.1,
                     encut=encut,
                     kpts=kpts,
                     gamma=True).write_input(cell)
            os.chdir('../')
            os.chdir('../')
            os.chdir('..')



                
        