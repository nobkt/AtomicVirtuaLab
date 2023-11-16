from AtomicVirtuaLab.build import interface_from_slab, slabgen
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.lammps import mk_nvt_input_uff
from AtomicVirtuaLab.tools import add_displacement, deform_cell
from AtomicVirtuaLab.io import rd_lammpsdata
from AtomicVirtuaLab.gpaw import mk_gpaw_pw_input_npt
from ase.io import read
from ase.visualize import view
from ase.build import make_supercell
import os
import sys

g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'

os.makedirs('interface_test',exist_ok=True)
os.chdir('interface_test')

P=[0]
T=[300,1000,2000,3000]
magElm="{'Fe':2.21}"

mpid = 2815

substrate_cell0 = read(g.cifdir+'/MoS2_mp'+str(mpid)+'.cif')
substrate_cell0.write('./bulk.cif')
os.system('atomsk bulk.cif -orthogonal-cell ortho.cfg')
substrate_cell0 = read('./ortho.cfg')
os.system('rm bulk.cif')
os.system('rm ortho.cfg')

film_cell0 = read(g.cifdir+'/Fe.cif')
film_cell0.write('./bulk.cif')
os.system('atomsk bulk.cif -orthogonal-cell ortho.cfg')
film_cell0 = read('./ortho.cfg')
os.system('rm bulk.cif')
os.system('rm ortho.cfg')

#view(substrate_cell0)
#view(film_cell0)

#substrate_slab_test = slabgen(substrate_cell0,0,0,1,1,1,1,7.5,7.5,fname='substrate',inv=False)
#film_slab_test = slabgen(film_cell0,1,0,0,1,1,1,7.5,7.5,fname='film',inv=False)
#view(substrate_slab_test)
#view(film_slab_test)
#sys.exit()

os.makedirs(str(mpid),exist_ok=True)
os.chdir(str(mpid))
interface_slab = interface_from_slab(substrate_cell0,0,0,0,1,3,2,2,film_cell0,0,1,0,0,3,3,3,in_plane_offset=(0.5,0.0),gap=1.6,vacuum_over_film=1.6)
print(len(interface_slab))
interface_slab.write('test.cif')
view(interface_slab)

for P0 in P:
    for T0 in T:
        cell_solid = interface_slab.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        #view(cell_solid)
        #sys.exit()
        cell_liquid = interface_slab.copy()
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #view(cell_liquid)
        #sys.exit()
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_gpaw_pw_input_npt(cell_solid,T0,P0,500,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=4.0,dftd3=True,berendsen=True,magElm=magElm)
        #mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',ecutwfc=ecutwfc0,ecutrho=ecutrho0,estep=9999,nstep=500,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=True,magElm=['Fe'])
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,10000.0,12345)
        os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        view(cell_liquid)
        os.chdir('../')
        os.system('rm -rf tmp')
        mk_gpaw_pw_input_npt(cell_liquid,T0,P0,500,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=4.0,dftd3=True,berendsen=True,magElm=magElm)
        #mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',ecutwfc=ecutwfc0,ecutrho=ecutrho0,estep=9999,nstep=500,ecut='manual',options={'vdw_corr':'dft-d3','dftd3_version':4},nspin=True,magElm=['Fe'])
        os.chdir('../')

os.chdir('../')