from AtomicVirtuaLab.espresso import mk_qe_input_scf, mk_qe_input_vcrelax
from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.io import read
from ase.build import make_supercell
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/inorgPigment'

os.makedirs('inorg_pigment',exist_ok=True)
os.chdir('inorg_pigment')

# opt U prame

Ulist = [5]

os.makedirs('test',exist_ok=True)
os.chdir('test')

os.makedirs('NBPM26',exist_ok=True)
os.chdir('NBPM26')

cell = rd_cif(g.cifs+'/NBPM26_final_9data_97-96-95-94-93-91-90-89-88.cif',primitive_cell=False)

kx = 4
ky = 4
kz = 4

os.makedirs('NM',exist_ok=True)
os.chdir('NM')
for U in Ulist:
    os.makedirs('U='+str(U),exist_ok=True)
    os.chdir('U='+str(U))
    #mk_qe_input_scf(cell,'pbe','paw',level='SSSP_precision',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='auto',tstress=True,nosym=True,options={},nspin=True)
    mk_qe_input_vcrelax(cell,'pbe','paw',level='SSSP_precision',estep=1000,nstep=1000,nosym=True,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='auto',options={},nspin=False)
    if U != 0:
        f = open('qe_vc-relax.pwi','a')
        f.write('HUBBARD (ortho-atomic)'+'\n')
        f.write('U Mn-3d '+str(float(U))+'\n')
        f.write('\n')
        f.close()
    os.chdir('../')
os.chdir('../')

os.makedirs('FM',exist_ok=True)
os.chdir('FM')
for U in Ulist:
    os.makedirs('U='+str(U),exist_ok=True)
    os.chdir('U='+str(U))
    #mk_qe_input_scf(cell,'pbe','paw',level='SSSP_precision',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='auto',tstress=True,nosym=True,options={},nspin=True)
    mk_qe_input_vcrelax(cell,'pbe','paw',level='SSSP_precision',estep=1000,nstep=1000,nosym=True,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=(kx,ky,kz),ecut='auto',options={},nspin=True)
    if U != 0:
        f = open('qe_vc-relax.pwi','a')
        f.write('HUBBARD (ortho-atomic)'+'\n')
        f.write('U Mn-3d '+str(float(U))+'\n')
        f.write('\n')
        f.close()
    os.chdir('../')