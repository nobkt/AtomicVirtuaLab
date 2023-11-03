from AtomicVirtuaLab.espresso import mk_qe_input_npt
from AtomicVirtuaLab.lammps import mk_nvt_input_uff
from AtomicVirtuaLab.io import rd_cif, rd_lammpsdata
from AtomicVirtuaLab.build import slabgen
from AtomicVirtuaLab.gpaw import mk_gpaw_pw_input_npt
from ase.build import make_supercell
from ase.visualize import view
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.tools import add_displacement, deform_cell, scale_cell
import os
import sys

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

#P = [0.0, 50.0, -50.0]
#T=[300,1000,2000,3000]


# train GPAW
#bulk
mpid = '2815'
os.makedirs('./MoS2/data/mp-'+str(mpid)+'/bulk',exist_ok=True)
os.chdir('./MoS2/data/mp-'+str(mpid)+'/bulk')

nx=3
ny=3
nz=2

cell = rd_cif(g.cifs+'/'+'MoS2_mp'+str(mpid)+'.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)
#view(cell)
#sys.exit()

P = [0.0,50000.0,-50000.0]
T=[300,1000,2000,3000,4000,5000]

for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        cell_liquid = cell.copy()
        cell_liquid = scale_cell(cell_liquid,0.8)
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        #cell_solid.write('cell.cif')
        mk_gpaw_pw_input_npt(cell_solid,T0,P0,500,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=4.0,dftd3=True,berendsen=True)
        #mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,200000.0,12345)
        os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        os.chdir('../')
        os.system('rm -rf tmp')
        #cell_liquid.write('cell.cif')
        mk_gpaw_pw_input_npt(cell_liquid,T0,P0,500,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=4.0,dftd3=True,berendsen=True)
        #mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
# train GPAW 終了



"""
# slab
mpid = 1434
nx001 = 3
ny001 = 3
nz001 = 1
nx100 = 5
ny100 = 1
nz100 = 4


cell = rd_cif(g.cifs+'/MoS2_mp1434.cif')

os.makedirs('./MoS2/data/mp-1434/slab',exist_ok=True)
os.chdir('./MoS2/data/mp-1434/slab')

os.makedirs('./001',exist_ok=True)
os.chdir('./001')
#slab = slabgen(cell,0,0,1,1,1,1,7.5,7.5)
#view(slab)
slab = slabgen(cell,0,0,1,nx001,ny001,nz001,7.5,7.5)
id = 2
#view(slab[id])

for P0 in P:
    for T0 in T:
        slab0 = slab[id].copy()
        slab0 = add_displacement(slab0,0.05)
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(slab0,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',cell_dofree='c',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
os.chdir('../')

os.makedirs('./100',exist_ok=True)
os.chdir('./100')
slab = slabgen(cell,1,0,0,1,1,1,7.5,7.5)
#view(slab)
#sys.exit()
slab = slabgen(cell,1,0,0,nx100,ny100,nz100,7.5,7.5)
id = 0
#view(slab[id])

for P0 in P:
    for T0 in T:
        slab0 = slab[id].copy()
        slab0 = add_displacement(slab0,0.05)
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(slab0,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',cell_dofree='c',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
os.chdir('../')
os.chdir('../')
"""

"""
mpid = 2815
nx001 = 3
ny001 = 3
nz001 = 2
nx100 = 5
ny100 = 2
nz100 = 4

cell = rd_cif(g.cifs+'/MoS2_mp2815.cif')

os.makedirs('./MoS2/data/mp-2815/slab',exist_ok=True)
os.chdir('./MoS2/data/mp-2815/slab')

os.makedirs('./001',exist_ok=True)
os.chdir('./001')
#slab = slabgen(cell,0,0,1,1,1,1,7.5,7.5)
#view(slab)
#exit()
slab = slabgen(cell,0,0,1,nx001,ny001,nz001,7.5,7.5)
id = 0
#view(slab[id])

for P0 in P:
    for T0 in T:
        slab0 = slab[id].copy()
        slab0 = add_displacement(slab0,0.05)
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(slab0,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',cell_dofree='c',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
os.chdir('../')

os.makedirs('./100',exist_ok=True)
os.chdir('./100')
#slab = slabgen(cell,1,0,0,1,1,1,7.5,7.5)
#view(slab)
#sys.exit()
slab = slabgen(cell,1,0,0,nx100,ny100,nz100,7.5,7.5)
id = 0
#view(slab[id])
#print(len(slab[0]))
#sys.exit()

for P0 in P:
    for T0 in T:
        slab0 = slab[id].copy()
        slab0 = add_displacement(slab0,0.05)
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(slab0,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',cell_dofree='c',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
os.chdir('../')
os.chdir('../')
"""

"""
#bulk
# mp2815
os.makedirs('./MoS2/data/mp-2815/bulk',exist_ok=True)
os.chdir('./MoS2/data/mp-2815/bulk')

nx=3
ny=3
nz=2

cell = rd_cif(g.cifs+'/'+'MoS2_mp2815.cif')
cell = make_supercell(cell,([nx,0,0],[0,ny,0],[0,0,nz]),wrap=True)
view(cell)
#print(len(cell))
#sys.exit()
P = [50.0, -50.0]
T=[300,1000,2000,3000]

for P0 in P:
    for T0 in T:
        cell_solid = cell.copy()
        cell_solid = add_displacement(cell_solid,0.05)
        cell_solid = deform_cell(cell_solid,0.05,0.5)
        cell_liquid = cell.copy()
        cell_liquid = scale_cell(cell_liquid,0.8)
        cell_liquid = add_displacement(cell_liquid,0.05)
        cell_liquid = deform_cell(cell_liquid,0.05,0.5)
        #
        # solid model
        #
        os.makedirs('./solid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./solid_T'+str(T0)+'_P'+str(P0))
        mk_qe_input_npt(cell_solid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
        #
        # liquid model
        #
        os.makedirs('./liquid_T'+str(T0)+'_P'+str(P0),exist_ok=True)
        os.chdir('./liquid_T'+str(T0)+'_P'+str(P0))
        os.makedirs('./tmp',exist_ok=True)
        os.chdir('./tmp')
        mk_nvt_input_uff(cell_liquid,0.0005,200,200,20000,200000.0,12345)
        os.system('lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
        cell_liquid = rd_lammpsdata(cell_liquid,'./result.data',True)
        os.chdir('../')
        os.system('rm -rf tmp')
        mk_qe_input_npt(cell_liquid,'pbe','paw',float(T0),100.0,float(P0),dt=4.0,level='high',estep=9999,nstep=500,ecut='auto',options={'vdw_corr':'dft-d3','dftd3_version':4})
        os.chdir('../')
"""

    
    
    
    
    

    
    
    
