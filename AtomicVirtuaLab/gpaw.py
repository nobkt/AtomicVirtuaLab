def mk_gpaw_pw_input_npt(cell0,temp,press,nstep,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=False,berendsen=False,berandsen_nstep=100,restart=False,magElm=None):
    f = open('gpw_npt.py','w')
    f.write('from ase.io import read'+'\n')
    f.write('from gpaw import GPAW, PW, Mixer, FermiDirac'+'\n')
    f.write('from ase.io import Trajectory'+'\n')
    f.write('from ase.md.npt import NPT'+'\n')
    f.write('from ase.md.velocitydistribution import MaxwellBoltzmannDistribution'+'\n')
    if berendsen:
        f.write('from ase.md.nptberendsen import NPTBerendsen'+'\n')
    f.write('import ase.units as units'+'\n')
    f.write('import math'+'\n')
    if dftd3:
        f.write('from dftd3.ase import DFTD3'+'\n')
    if restart:
        f.write("cell = read('opt.traj')"+'\n')
    else:
        cell0.write('cell.cif')
    if not restart:
        f.write('cell = read("cell.cif")'+'\n')
    if magElm is not None:
        f.write('magmom=[0.0]*int(len(cell))'+'\n')
        f.write('magElm = '+str(magElm)+'\n')
        f.write('for elm in magElm:'+'\n')
        f.write('    for atom in cell:'+'\n')
        f.write('        if atom.symbol == elm:'+'\n')
        f.write('            magmom.append(float(magElm[elm]))'+'\n')
        f.write('cell.set_initial_magnetic_moments(magmom)'+'\n')
    f.write("dft = GPAW(mode=PW("+str(ecut)+"),xc='"+str(xc)+"',kpts="+str(kpts)+",maxiter="+str(maxiter)+",occupations=FermiDirac(width=0.05),mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0),txt='npt_gs.txt')"+'\n')
    if dftd3:
        f.write("d3 = DFTD3(method='"+str(xc)+"',damping='d3bj')"+'\n')
        f.write('cell.calc = d3.add_calculator(dft)'+'\n')
    else:
        f.write('cell.calc = dft'+'\n')
    f.write('cell.get_potential_energy()'+'\n')
    f.write('MaxwellBoltzmannDistribution(cell, temperature_K='+str(temp)+')'+'\n')
    if berendsen:
        f.write('dyn = NPTBerendsen(atoms=cell,timestep='+str(dt)+'*units.fs,temperature_K='+str(temp)+',pressure_au='+str(press)+'*units.bar,taut=100*units.fs,taup=1000.0*units.fs,compressibility_au=4.57e-5/units.bar,logfile="log_berendsen_npt",trajectory="berendsen_npt.traj",loginterval=1)'+'\n')
        f.write('dyn.run('+str(berandsen_nstep)+')'+'\n')
    f.write('dyn = NPT(atoms=cell,timestep='+str(dt)+'*units.fs,temperature_K='+str(temp)+',externalstress='+str(press)+'*units.bar,ttime='+str(25.0)+'*units.fs,pfactor='+str(75.0*75.0*100.0)+'*units.GPa*(units.fs**2),logfile="log_pr_npt",trajectory="pr_npt.traj",loginterval=1)'+'\n')
    f.write('dyn.run('+str(nstep)+')'+'\n')
    f.close()

def mk_gpaw_lcao_input_npt(cell0,temp,press,nstep,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dt=0.5,dftd3=False,berendsen=False,berandsen_nstep=100,restart=False):
    f = open('gpw_npt.py','w')
    f.write('from ase.io import read'+'\n')
    f.write('from gpaw import GPAW, PW, Mixer, FermiDirac'+'\n')
    f.write('from ase.io import Trajectory'+'\n')
    f.write('from ase.md.npt import NPT'+'\n')
    f.write('from ase.md.velocitydistribution import MaxwellBoltzmannDistribution'+'\n')
    if berendsen:
        f.write('from ase.md.nptberendsen import NPTBerendsen'+'\n')
    f.write('import ase.units as units'+'\n')
    f.write('import math'+'\n')
    if dftd3:
        f.write('from dftd3.ase import DFTD3'+'\n')
    if restart:
        f.write("cell = read('opt.traj')"+'\n')
    else:
        cell0.write('cell.cif')
    if not restart:
        f.write('cell = read("cell.cif")'+'\n')
    f.write("dft = GPAW(mode='lcao',basis='"+str(basis)+"',xc='"+str(xc)+"',kpts="+str(kpts)+",maxiter="+str(maxiter)+",occupations=FermiDirac(width=0.05),mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0),txt='npt_gs.txt')"+'\n')
    if dftd3:
        f.write("d3 = DFTD3(method='"+str(xc)+"',damping='d3bj')"+'\n')
        f.write('cell.calc = d3.add_calculator(dft)'+'\n')
    else:
        f.write('cell.calc = dft'+'\n')
    f.write('cell.get_potential_energy()'+'\n')
    f.write('MaxwellBoltzmannDistribution(cell, temperature_K='+str(temp)+')'+'\n')
    if berendsen:
        f.write('dyn = NPTBerendsen(atoms=cell,timestep='+str(dt)+'*units.fs,temperature_K='+str(temp)+',pressure_au='+str(press)+'*units.bar,taut=100*units.fs,taup=1000.0*units.fs,compressibility_au=4.57e-5/units.bar,logfile="log_berendsen_npt",trajectory="berendsen_npt.traj",loginterval=1)'+'\n')
        f.write('dyn.run('+str(berandsen_nstep)+')'+'\n')
    f.write('dyn = NPT(atoms=cell,timestep='+str(dt)+'*units.fs,temperature_K='+str(temp)+',externalstress='+str(press)+'*units.bar,ttime='+str(25.0)+'*units.fs,pfactor='+str(75.0*75.0*100.0)+'*units.GPa*(units.fs**2),logfile="log_pr_npt",trajectory="pr_npt.traj",loginterval=1)'+'\n')
    f.write('dyn.run('+str(nstep)+')'+'\n')
    f.close()

def mk_gpaw_pw_input_optimize(cell0,xc='PBE',ecut=400,kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=False):
    cell0.write('cell.cif')
    f = open('gpw_opt.py','w')
    f.write('from ase.io import read'+'\n')
    f.write('from gpaw import GPAW, PW, Mixer, FermiDirac'+'\n')
    f.write('from ase.io import Trajectory'+'\n')
    f.write('from ase.optimize import QuasiNewton'+'\n')
    f.write('import ase.units as units'+'\n')
    f.write('import math'+'\n')
    if dftd3:
        f.write('from dftd3.ase import DFTD3'+'\n')
    f.write('cell = read("cell.cif")'+'\n')
    f.write("dft = GPAW(mode=PW("+str(ecut)+"),xc='"+str(xc)+"',kpts="+str(kpts)+",maxiter="+str(maxiter)+",occupations=FermiDirac(width=0.05),mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0),txt='opt_gs.txt')"+'\n')
    if dftd3:
        f.write("d3 = DFTD3(method='"+str(xc)+"',damping='d3bj')"+'\n')
        f.write('cell.calc = d3.add_calculator(dft)'+'\n')
    else:
        f.write('cell.calc = dft'+'\n')
    f.write('cell.get_potential_energy()'+'\n')
    f.write("dyn = QuasiNewton(cell,trajectory='opt.traj')"+'\n')
    f.write('dyn.run(fmax=0.05)'+'\n')
    f.close()

def mk_gpaw_lcao_input_optimize(cell0,xc='PBE',basis='dzp',kpts={'size':(1,1,1),'gamma':True},maxiter=2000,dftd3=False):
    cell0.write('cell.cif')
    f = open('gpw_opt.py','w')
    f.write('from ase.io import read'+'\n')
    f.write('from gpaw import GPAW, PW, Mixer, FermiDirac'+'\n')
    f.write('from ase.io import Trajectory'+'\n')
    f.write('from ase.optimize import QuasiNewton'+'\n')
    f.write('import ase.units as units'+'\n')
    f.write('import math'+'\n')
    if dftd3:
        f.write('from dftd3.ase import DFTD3'+'\n')
    f.write('cell = read("cell.cif")'+'\n')
    f.write("dft = GPAW(mode='lcao',basis='"+str(basis)+"',xc='"+str(xc)+"',kpts="+str(kpts)+",maxiter="+str(maxiter)+",occupations=FermiDirac(width=0.05),mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0),txt='opt_gs.txt')"+'\n')
    if dftd3:
        f.write("d3 = DFTD3(method='"+str(xc)+"',damping='d3bj')"+'\n')
        f.write('cell.calc = d3.add_calculator(dft)'+'\n')
    else:
        f.write('cell.calc = dft'+'\n')
    f.write('cell.get_potential_energy()'+'\n')
    f.write("dyn = QuasiNewton(cell,trajectory='opt.traj')"+'\n')
    f.write('dyn.run(fmax=0.5)'+'\n')
    f.close()
        
    