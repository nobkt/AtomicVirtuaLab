def mk_nvt_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,md_step,T,seed,restart):
    f = open('input_nvt.lmp','w')
    f.write('# ------------------------------- Initialization Section -------------------'+'\n')
    f.write('\n')
    f.write('include "system.in.init"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Atom Definition Section -------------------'+'\n')
    f.write('\n')
    if restart:
        f.write('read_data "restart.data"'+'\n')
    else:
        f.write('read_data "system.data"'+'\n')
        f.write('\n')
        f.write('include "system.in.charges"'+'\n')
        f.write('\n')
        f.write('# ------------------------------- Settings Section --------------------------'+'\n')
        f.write('\n')
        f.write('include "system.in.settings"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Run Section -------------------------------'+'\n')
    f.write('\n')
    if not restart:
        if minimize:
            f.write('# -- minimization protocol --'+'\n')
            f.write('\n')
            f.write('minimize 1.0e-4 1.0e-6 100000 400000'+'\n')
            f.write('\n')
    f.write('# -- simulation protocol --'+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('print "--- Run a simulation using a Nose-Hoover Thermostat ---"'+'\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    if not restart:
        f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnvt all nvt temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data system_after_nvt.data'+'\n')
    f.close()
    
def mk_npt_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,md_step,T,P,dof_type,seed,restart,qeq=False):
    f = open('input_npt.lmp','w')
    f.write('# ------------------------------- Initialization Section -------------------'+'\n')
    f.write('\n')
    f.write('include "system.in.init"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Atom Definition Section -------------------'+'\n')
    f.write('\n')
    if restart:
        f.write('read_data "restart.data"'+'\n')
    else:
        f.write('read_data "system.data"'+'\n')
        f.write('\n')
        if qeq:
            f.write('set atom 1 charge 0.0001'+'\n')
            f.write('set atom 2 charge -0.0001'+'\n')
        else:
            f.write('include "system.in.charges"'+'\n')
        f.write('\n')
        f.write('# ------------------------------- Settings Section --------------------------'+'\n')
        f.write('\n')
        f.write('include "system.in.settings"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Run Section -------------------------------'+'\n')
    f.write('\n')
    if qeq:
        mk_qeqfile(symbols,unit='real')
        f.write('fix fxqeq all qeq/point 1 10 1.0e-6 100 my_qeq'+'\n')
    if not restart:
        if minimize:
            f.write('# -- minimization protocol --'+'\n')
            f.write('\n')
            f.write('minimize 1.0e-4 1.0e-6 100000 400000'+'\n')
            f.write('\n')
    f.write('# -- simulation protocol --'+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('print "--- Run a simulation using a Nose-Hoover Thermostat/Barostat ---"'+'\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    if not restart:
        f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnpt all npt temp '+str(T)+' '+str(T)+' $(dt*100.0) '+str(dof_type)+' '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data system_after_npt.data'+'\n')
    f.close()

def mk_npt_compress_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,ncomp,comp_step,Tcomp,Pcomp,hT_step,hT,md_step,T,P,dof_type,seed,restart,qeq=False,rdfpairs=None):
    f = open('input_npt.lmp','w')
    f.write('# ------------------------------- Initialization Section -------------------'+'\n')
    f.write('\n')
    f.write('include "system.in.init"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Atom Definition Section -------------------'+'\n')
    f.write('\n')
    if restart:
        f.write('read_data "restart.data"'+'\n')
    else:
        f.write('read_data "system.data"'+'\n')
        f.write('\n')
        if qeq:
            f.write('set atom 1 charge 0.0001'+'\n')
            f.write('set atom 2 charge -0.0001'+'\n')
        else:
            f.write('include "system.in.charges"'+'\n')
        f.write('\n')
        f.write('# ------------------------------- Settings Section --------------------------'+'\n')
        f.write('\n')
        f.write('include "system.in.settings"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Run Section -------------------------------'+'\n')
    f.write('\n')
    if qeq:
        mk_qeqfile(symbols,unit='real')
        f.write('fix fxqeq all qeq/point 1 10 1.0e-6 100 my_qeq'+'\n')
    f.write('fix fxmom all momentum 1 linear 1 1 1'+'\n')
    if not restart:
        if minimize:
            f.write('# -- minimization protocol --'+'\n')
            f.write('\n')
            f.write('minimize 1.0e-4 1.0e-6 100000 400000'+'\n')
            f.write('\n')
    f.write('# -- simulation protocol --'+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('print "--- Run a simulation using a Nose-Hoover Thermostat/Barostat ---"'+'\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_comp.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    if not restart:
        f.write('velocity all create '+str(Tcomp)+' '+str(seed)+'\n')
    f.write('fix fxnpt all npt temp '+str(Tcomp)+' '+str(Tcomp)+' $(dt*100.0) '+str(dof_type)+' '+str(Pcomp)+' '+str(Pcomp)+' $(dt*1000.0)'+'\n')
    f.write('\n')
    f.write('run '+str(comp_step)+'\n')
    f.write('unfix fxnpt'+'\n')
    f.write('velocity all scale '+str(Tcomp)+'\n')
    f.write('fix fxnvt all nvt temp '+str(hT)+' '+str(hT)+' $(dt*100.0)'+'\n')
    f.write('\n')
    f.write('run '+str(hT_step)+'\n')
    f.write('unfix fxnvt'+'\n')
    #
    f.write('label loopa'+'\n')
    f.write('variable a loop '+str(ncomp-1)+'\n')
    f.write('fix fxnpt all npt temp '+str(Tcomp)+' '+str(Tcomp)+' $(dt*100.0) '+str(dof_type)+' '+str(Pcomp)+' '+str(Pcomp)+' $(dt*1000.0)'+'\n')
    f.write('\n')
    f.write('run '+str(comp_step)+'\n')
    f.write('unfix fxnpt'+'\n')
    f.write('velocity all scale '+str(Tcomp)+'\n')
    f.write('fix fxnvt all nvt temp '+str(hT)+' '+str(hT)+' $(dt*100.0)'+'\n')
    f.write('\n')
    f.write('run '+str(hT_step)+'\n')
    f.write('unfix fxnvt'+'\n')
    f.write('next a'+'\n')
    f.write('jump SELF loopa'+'\n')
    f.write('undump dmp'+'\n')
    #
    f.write('reset_timestep 0'+'\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    if rdfpairs is not None:
        f.write('compute myRDF all rdf 1000')
        for rdfpair in rdfpairs:
            f.write(' '+str(rdfpair[0])+' '+str(rdfpair[1]))
        f.write('\n')
        f.write('fix trdf all ave/time 1 '+str(int(md_step/2))+' '+str(md_step)+' c_myRDF[*] file tmp.rdf mode vector'+'\n')
    f.write('fix th all ave/time 1 '+str(int(md_step/2))+' '+str(md_step)+' v_mytemp v_myenthalpy v_mydensity v_myvol v_mycella v_mycellb v_mycellc v_mycellalpha v_mycellbeta v_mycellgamma file th.profile'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('velocity all scale '+str(T)+'\n')
    f.write('fix fxnpt all npt temp '+str(T)+' '+str(T)+' $(dt*100.0) '+str(dof_type)+' '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data system_after_npt.data'+'\n')
    f.close()

def mk_nvtdeform_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,md_step,T,erate,tensile_step,seed,restart):
    f = open('input_nvtdeform.lmp','w')
    f.write('# ------------------------------- Initialization Section -------------------'+'\n')
    f.write('\n')
    f.write('include "system.in.init"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Atom Definition Section -------------------'+'\n')
    f.write('\n')
    if restart:
        f.write('read_data "restart.data"'+'\n')
    else:
        f.write('read_data "system.data"'+'\n')
        f.write('\n')
        f.write('include "system.in.charges"'+'\n')
        f.write('\n')
        f.write('# ------------------------------- Settings Section --------------------------'+'\n')
        f.write('\n')
        f.write('include "system.in.settings"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Run Section -------------------------------'+'\n')
    f.write('\n')
    if not restart:
        if minimize:
            f.write('# -- minimization protocol --'+'\n')
            f.write('\n')
            f.write('minimize 1.0e-4 1.0e-6 100000 400000'+'\n')
            f.write('\n')
    f.write('# -- simulation protocol --'+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('print "--- Run a simulation using a Nose-Hoover Thermostat + Defrom---"'+'\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvtdeform.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('compute stressAtom all stress/atom NULL'+'\n')
    f.write('compute myStress all reduce sum c_stressAtom[3]'+'\n')
    f.write('compute pp all pressure NULL virial'+'\n')
    f.write('variable myStress equal c_myStress/vol*0.000101325'+'\n')
    f.write('variable pzz equal -1.0*c_pp[3]*0.000101325'+'\n')
    f.write('variable tmp equal "lz"'+'\n')
    f.write('variable L0 equal ${tmp}'+'\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma v_myStress v_pzz'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    if not restart:
        f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnvt all nvt temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('fix fxdeform all deform 1 z erate '+str(erate)+' units box remap x'+'\n')
    f.write('variable strain equal "(lz-v_L0)/v_L0"'+'\n')
    f.write('variable p1 equal "v_strain"'+'\n')
    f.write('variable p2 equal "v_myStress"'+'\n')
    f.write('variable p3 equal "v_pzz"'+'\n')
    f.write('fix fxprint all print '+str(tensile_step)+' "${p1} ${p2} ${p3}" file tensile.txt screen no'+'\n')
    f.write('\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data system_after_nvtdeform.data'+'\n')
    f.close()

def mk_nptdeform_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,md_step,T,P,deform_style,erate,tensile_step,seed,restart):
    f = open('input_nvtdeform.lmp','w')
    f.write('# ------------------------------- Initialization Section -------------------'+'\n')
    f.write('\n')
    f.write('include "system.in.init"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Atom Definition Section -------------------'+'\n')
    f.write('\n')
    if restart:
        f.write('read_data "restart.data"'+'\n')
    else:
        f.write('read_data "system.data"'+'\n')
        f.write('\n')
        f.write('include "system.in.charges"'+'\n')
        f.write('\n')
        f.write('# ------------------------------- Settings Section --------------------------'+'\n')
        f.write('\n')
        f.write('include "system.in.settings"'+'\n')
    f.write('\n')
    f.write('# ------------------------------- Run Section -------------------------------'+'\n')
    f.write('\n')
    if not restart:
        if minimize:
            f.write('# -- minimization protocol --'+'\n')
            f.write('\n')
            f.write('minimize 1.0e-4 1.0e-6 100000 400000'+'\n')
            f.write('\n')
    f.write('# -- simulation protocol --'+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('print "--- Run a simulation using a Nose-Hoover Thermostat/Barostat + Defrom---"'+'\n')
    f.write('print "----------------------------------------------------------------"'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nptdeform.lammpstrj id mol type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('compute stressAtom all stress/atom NULL'+'\n')
    f.write('compute myStress all reduce sum c_stressAtom[3]'+'\n')
    f.write('compute pp all pressure NULL virial'+'\n')
    f.write('variable myStress equal c_myStress/vol*0.000101325'+'\n')
    f.write('variable pzz equal -1.0*c_pp[3]*0.000101325'+'\n')
    f.write('variable tmp equal "lz"'+'\n')
    f.write('variable L0 equal ${tmp}'+'\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma v_myStress v_pzz'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    if not restart:
        f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    if deform_style == 'semiisotropic':
        f.write('fix fxnpt all npt temp '+str(T)+' '+str(T)+' $(100.0*dt) '+'x'+' '+str(P)+' '+str(P)+' $(1000.0*dt) '+'y'+' '+str(P)+' '+str(P)+' $(1000.0*dt) '+'couple xy'+'\n')
    elif deform_style == 'anisotorpic':
        f.write('fix fxnpt all npt temp '+str(T)+' '+str(T)+' $(100.0*dt) '+'x'+' '+str(P)+' '+str(P)+' $(1000.0*dt) '+'y'+' '+str(P)+' '+str(P)+' $(1000.0*dt)'+'\n')
    f.write('fix fxdeform all deform 1 z erate '+str(erate)+' units box remap x'+'\n')
    f.write('variable strain equal "(lz-v_L0)/v_L0"'+'\n')
    f.write('variable p1 equal "v_strain"'+'\n')
    f.write('variable p2 equal "v_myStress"'+'\n')
    f.write('variable p3 equal "v_pzz"'+'\n')
    f.write('fix fxprint all print '+str(tensile_step)+' "${p1} ${p2} ${p3}" file tensile.txt screen no'+'\n')
    f.write('\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data system_after_nptdeform.data'+'\n')
    f.close()

def mk_nvt_input_uff(cell,dt,dmp_step,thermo_step,md_step,T,seed):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    mk_qeqfile(symbols)
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style charge'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('set atom 1 charge 0.0001'+'\n')
    f.write('set atom 2 charge -0.0001'+'\n')
    f.write('pair_style lj/cut/coul/long 10.0 10.0'+'\n')
    f.write('kspace_style ewald 0.000001'+'\n')
    i = 1
    for symbol in symbols:
        sigma, epsilon = set_uff_lj(symbol)
        f.write('pair_coeff '+str(i)+' '+str(i)+' '+str(epsilon)+' '+str(sigma)+'\n')
        i = i + 1
    f.write('pair_modify mix geometric'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('fix fxqeq all qeq/point 1 10 1.0e-6 100 my_qeq'+'\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnvt all nvt temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data result.data'+'\n')
    f.close()

def mk_mimize_input_uff_adsorp(cell,nstep,dmp_step,thermo_step):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,False,mol=True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p f'+'\n')
    f.write('atom_style full'+'\n')
    f.write('read_data lammps.data'+'\n')
    f.write('group mol molecule 2'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style lj/cut 10.0'+'\n')
    i = 1
    for symbol in symbols:
        sigma, epsilon = set_uff_lj(symbol)
        f.write('pair_coeff '+str(i)+' '+str(i)+' '+str(epsilon)+' '+str(sigma)+'\n')
        i = i + 1
    f.write('pair_modify mix geometric'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('velocity mol create 5.0 12345'+'\n')
    f.write('fix fxnvt mol rigid/nvt molecule temp 5.0 5.0 $(dt*100.0)'+'\n')
    f.write('run '+str(nstep)+'\n')
    f.write('write_data result.data'+'\n')
    f.write('variable mype equal pe'+'\n')
    f.write('run 0'+'\n')
    f.write('print "Potential Energy: ${mype}"'+'\n')
    f.close()

def mk_mimize_input_uff_adsorp_onMoS2(cell,nstep,dmp_step,thermo_step):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,False,mol=True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p f'+'\n')
    f.write('atom_style full'+'\n')
    f.write('read_data lammps.data'+'\n')
    f.write('group mol molecule 4'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style lj/cut 10.0'+'\n')
    i = 1
    for symbol in symbols:
        sigma, epsilon = set_uff_lj(symbol)
        f.write('pair_coeff '+str(i)+' '+str(i)+' '+str(epsilon)+' '+str(sigma)+'\n')
        i = i + 1
    f.write('pair_modify mix geometric'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('velocity mol create 5.0 12345'+'\n')
    f.write('fix fxnvt mol rigid/nvt molecule temp 5.0 5.0 $(dt*100.0)'+'\n')
    f.write('run '+str(nstep)+'\n')
    f.write('write_data result.data'+'\n')
    f.write('variable mype equal pe'+'\n')
    f.write('run 0'+'\n')
    f.write('print "Potential Energy: ${mype}"'+'\n')
    f.close()

def mk_mimize_input_dp_adsorp_onMoS2(cell,nstep,dmp_step,thermo_step):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,False,mol=True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p f'+'\n')
    f.write('atom_style full'+'\n')
    f.write('read_data lammps.data'+'\n')
    f.write('group mol molecule 4'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('velocity mol create 5.0 12345'+'\n')
    f.write('fix fxnvt mol nvt temp 5.0 5.0 $(dt*100.0)'+'\n')
    f.write('run '+str(nstep)+'\n')
    f.write('minimize 1.0e-10 1.0e-10 10000 100000'+'\n')
    f.write('write_data result.data'+'\n')
    f.write('variable mype equal pe'+'\n')
    f.write('run 0'+'\n')
    f.write('print "Potential Energy: ${mype}"'+'\n')
    f.close()

def mk_nvt_input_uff_rigid(cell,dt,dmp_step,thermo_step,md_step,T,seed):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    mk_qeqfile(symbols)
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style full'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('set atom 1 charge 0.0001'+'\n')
    f.write('set atom 2 charge -0.0001'+'\n')
    f.write('pair_style lj/cut/coul/long 10.0 10.0'+'\n')
    f.write('kspace_style ewald 0.000001'+'\n')
    i = 1
    for symbol in symbols:
        sigma, epsilon = set_uff_lj(symbol)
        f.write('pair_coeff '+str(i)+' '+str(i)+' '+str(epsilon)+' '+str(sigma)+'\n')
        i = i + 1
    f.write('pair_modify mix geometric'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('fix fxqeq all qeq/point 1 10 1.0e-6 100 my_qeq'+'\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnvt all rigid/nvt molecule temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data result.data'+'\n')
    f.close()

def mk_nvt_input_uff_rigid_scale(cell,dt,dmp_step,thermo_step,lat,df_step,md_step,T,seed):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,False,mol=True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    mk_qeqfile(symbols)
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style full'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('set atom 1 charge 0.0001'+'\n')
    f.write('set atom 2 charge -0.0001'+'\n')
    f.write('pair_style lj/cut/coul/long 10.0 10.0'+'\n')
    f.write('kspace_style ewald 0.000001'+'\n')
    i = 1
    for symbol in symbols:
        sigma, epsilon = set_uff_lj(symbol)
        f.write('pair_coeff '+str(i)+' '+str(i)+' '+str(epsilon)+' '+str(sigma)+'\n')
        i = i + 1
    f.write('pair_modify mix geometric'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('fix fxmom all momentum 1 linear 1 1 1'+'\n')
    f.write('fix fxqeq all qeq/point 1 10 1.0e-6 100 my_qeq'+'\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(T)+' '+str(seed)+'\n')
    f.write('fix fxnvt all rigid/nvt molecule temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('fix fdeform all deform 1 x final 0.0 '+str(lat)+' y final 0.0 '+str(lat)+' z final 0.0 '+str(lat)+' units box'+'\n')
    f.write('run '+str(df_step)+'\n')
    f.write('unfix fdeform'+'\n')
    f.write('run '+str(md_step)+'\n')
    f.write('write_data result.data'+'\n')
    f.close()

def mk_nvt_input_deepmd_friction(cell,dt,dmp_step,thermo_step,eq_step,md_step,nave,T,seed,nlayer=0,vdirec='y'):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    import random
    mk_lammpsdata(cell,False)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p f'+'\n')
    f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    f.write('region upper block INF INF INF INF ${z1} INF'+'\n')
    f.write('region lower block INF INF INF INF INF ${z0}'+'\n')
    f.write('group Mo type 1'+'\n')
    f.write('group upper region upper'+'\n')
    f.write('group lower region lower'+'\n')
    f.write('group upperMo intersect upper Mo'+'\n')
    f.write('group lowerMo intersect lower Mo'+'\n')
    
    if nlayer != 0:
        for nl in range(nlayer-2):
            if nl+1 == 13:
                break
            if nl == 0:
                f.write('region layer'+str(nl+1)+' block INF INF INF INF ${zm'+str(nl+1)+'} ${z1}'+'\n')
                f.write('group layer'+str(nl+1)+' region layer'+str(nl+1)+'\n')
                f.write('group layer'+str(nl+1)+'Mo intersect layer'+str(nl+1)+' Mo'+'\n')
            elif nl == nlayer-2-1:
                f.write('region layer'+str(nl+1)+' block INF INF INF INF ${z0} ${zm'+str(nl)+'}'+'\n')
                f.write('group layer'+str(nl+1)+' region layer'+str(nl+1)+'\n')
                f.write('group layer'+str(nl+1)+'Mo intersect layer'+str(nl+1)+' Mo'+'\n')
            else:
                f.write('region layer'+str(nl+1)+' block INF INF INF INF ${zm'+str(nl+1)+'} ${zm'+str(nl)+'}'+'\n')
                f.write('group layer'+str(nl+1)+' region layer'+str(nl+1)+'\n')
                f.write('group layer'+str(nl+1)+'Mo intersect layer'+str(nl+1)+' Mo'+'\n')
    
    f.write('group middle subtract all upper lower'+'\n')
    
    f.write('compute mype all pe/atom'+'\n')
    f.write('compute myke all ke/atom'+'\n')
    #f.write('compute sig all centroid/stress/atom NULL virial'+'\n')
    #if nlayer != 0:
    #    for nl in range(nlayer-2):
    #        if nl == 0:
    #            f.write('compute uptol'+str(nl+1)+' upper group/group layer'+str(nl+1)+'\n')
    #        elif nl == nlayer-2-1:
    #            f.write('compute l'+str(nl+1)+'tolow layer'+str(nl+1)+' group/group lower'+'\n')
    #        else:
    #            f.write('compute l'+str(nl)+'tol'+str(nl+1)+' layer'+str(nl)+' group/group layer'+str(nl+1)+'\n')
    
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    #f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt_load.lammpstrj id type element x y z vx vy vz fx fy fz c_sig[1] c_sig[2] c_sig[3] c_sig[4] c_sig[5] c_sig[6] c_sig[7] c_sig[8] c_sig[9] ix iy iz'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt_load.lammpstrj id type element x y z vx vy vz fx fy fz c_mype c_myke ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('fix rigid1 lower rigid single torque * off off off'+'\n')
    f.write('fix rigid2 upper rigid single torque * off off off'+'\n')
    f.write('variable S equal lx*ly'+'\n')
    f.write('variable P equal v_P0*1.0e-6/1.602176634'+'\n')
    f.write('variable N equal count(upper)'+'\n')
    f.write('variable Pz0 equal v_P*v_S'+'\n')
    f.write('variable Pz equal v_P*v_S/v_N'+'\n')
    f.write('fix freeze1 lower setforce 0.0 0.0 0.0'+'\n')
    f.write('fix freeze2 upper setforce 0.0 0.0 NULL'+'\n')
    f.write('fix addf upper addforce 0 0 v_Pz'+'\n')
    f.write('velocity middle create '+str(T)+' '+str(seed)+'\n')
    if vdirec == 'y':
        f.write('variable Fric equal f_freeze2[2]/v_Pz0'+'\n')
    elif vdirec == 'x':
        f.write('variable Fric equal f_freeze2[1]/v_Pz0'+'\n')

    #f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma f_freeze2[1] f_freeze2[2] f_freeze2[3] v_Pz0 v_Fric'+'\n')
    #f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('fix fxnvt middle nvt temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('run '+str(eq_step)+'\n')
    f.write('unfix freeze1'+'\n')
    f.write('unfix freeze2'+'\n')
    f.write('unfix rigid1'+'\n')
    f.write('unfix rigid2'+'\n')
    f.write('unfix fxnvt'+'\n')
    f.write('unfix addf'+'\n')
    f.write('undump dmp'+'\n')
    f.write('\n')

    f.write('dump dmp all custom '+str(dmp_step)+' traj_nvt_fric.lammpstrj id type element x y z vx vy vz fx fy fz c_mype c_myke ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    
    f.write('reset_timestep 0'+'\n')
    
    if nlayer != 0:
        f.write('compute uppermsd upper msd average yes'+'\n')
        f.write('compute lowermsd lower msd average yes'+'\n')
        f.write('compute upperpeperatom upper pe/atom'+'\n')
        f.write('compute lowerpeperatom lower pe/atom'+'\n')
        f.write('compute upperpe upper reduce sum c_upperpeperatom'+'\n')
        f.write('compute lowerpe lower reduce sum c_lowerpeperatom'+'\n')
        #f.write('compute upperpropperatom upper property/atom vx vy vz fx fy fz'+'\n')
        #f.write('compute lowerpropperatom lower property/atom vx vy vz fx fy fz'+'\n')
        #f.write('compute upperprop upper reduce ave c_upperpropperatom[*]'+'\n')
        #f.write('compute lowerprop lower reduce ave c_lowerpropperatom[*]'+'\n')
        f.write('compute upperMopropperatom upper property/atom z'+'\n')
        f.write('compute lowerMopropperatom lower property/atom z'+'\n')
        f.write('compute upperMomaxz upperMo reduce max c_upperMopropperatom'+'\n')
        f.write('compute lowerMomaxz lowerMo reduce max c_lowerMopropperatom'+'\n')
        f.write('compute upperMominz upperMo reduce min c_upperMopropperatom'+'\n')
        f.write('compute lowerMominz lowerMo reduce min c_lowerMopropperatom'+'\n')

    if nlayer != 0:
        for nl in range(nlayer-2):
            if nl+1 == 13:
                break
            if nl == 0:
                f.write('compute layer'+str(nl+1)+'msd layer'+str(nl+1)+' msd average yes'+'\n')
                f.write('compute layer'+str(nl+1)+'peperatom layer'+str(nl+1)+' pe/atom'+'\n')
                f.write('compute layer'+str(nl+1)+'pe layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'peperatom'+'\n')
                #f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom vx vy vz fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'prop layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'propperatom[*]'+'\n')
                f.write('compute layerMo'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom z'+'\n')
                f.write('compute layerMo'+str(nl+1)+'maxz layer'+str(nl+1)+' reduce max c_layerMo'+str(nl+1)+'propperatom'+'\n')
                f.write('compute layerMo'+str(nl+1)+'minz layer'+str(nl+1)+' reduce min c_layerMo'+str(nl+1)+'propperatom'+'\n')
            elif nl == nlayer-2-1:
                f.write('compute layer'+str(nl+1)+'msd layer'+str(nl+1)+' msd average yes'+'\n')
                f.write('compute layer'+str(nl+1)+'peperatom layer'+str(nl+1)+' pe/atom'+'\n')
                f.write('compute layer'+str(nl+1)+'pe layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'peperatom'+'\n')
                #f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom vx vy vz fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'prop layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'propperatom[*]'+'\n')
                f.write('compute layerMo'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom z'+'\n')
                f.write('compute layerMo'+str(nl+1)+'maxz layer'+str(nl+1)+' reduce max c_layerMo'+str(nl+1)+'propperatom'+'\n')
                f.write('compute layerMo'+str(nl+1)+'minz layer'+str(nl+1)+' reduce min c_layerMo'+str(nl+1)+'propperatom'+'\n')
            else:
                f.write('compute layer'+str(nl+1)+'msd layer'+str(nl+1)+' msd average yes'+'\n')
                f.write('compute layer'+str(nl+1)+'peperatom layer'+str(nl+1)+' pe/atom'+'\n')
                f.write('compute layer'+str(nl+1)+'pe layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'peperatom'+'\n')
                #f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom vx vy vz fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom fx fy fz'+'\n')
                f.write('compute layer'+str(nl+1)+'prop layer'+str(nl+1)+' reduce sum c_layer'+str(nl+1)+'propperatom[*]'+'\n')
                f.write('compute layerMo'+str(nl+1)+'propperatom layer'+str(nl+1)+' property/atom z'+'\n')
                f.write('compute layerMo'+str(nl+1)+'maxz layer'+str(nl+1)+' reduce max c_layerMo'+str(nl+1)+'propperatom'+'\n')
                f.write('compute layerMo'+str(nl+1)+'minz layer'+str(nl+1)+' reduce min c_layerMo'+str(nl+1)+'propperatom'+'\n')

    f.write('fix rigid1 lower rigid single torque * off off off'+'\n')
    f.write('fix rigid2 upper rigid single torque * off off off'+'\n')
    f.write('fix addf upper addforce 0 0 v_Pz'+'\n')
    f.write('variable vp equal ${vfric}'+'\n')
    f.write('variable vm equal -1.0*${vfric}'+'\n')
    if vdirec == 'y':
        f.write('fix freeze1 lower setforce NULL 0.0 0.0'+'\n')
        f.write('fix freeze2 upper setforce NULL 0.0 NULL'+'\n')
        f.write('velocity lower set NULL v_vm 0.0 units box'+'\n')
        f.write('velocity upper set NULL v_vp NULL units box'+'\n')
    elif vdirec == 'x':
        f.write('fix freeze1 lower setforce 0.0 NULL 0.0'+'\n')
        f.write('fix freeze2 upper setforce 0.0 NULL NULL'+'\n')
        f.write('velocity lower set v_vm NULL 0.0 units box'+'\n')
        f.write('velocity upper set v_vp NULL NULL units box'+'\n')
    f.write('fix th all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' f_freeze2[1] f_freeze2[2] f_freeze2[3] v_Pz0 v_Fric file friction.profile'+'\n')

    if nlayer != 0:
        for nl in range(nlayer-2):
            if nl+1 == 13:
                break
            if nl == 0:
                if vdirec == 'y':
                    f.write('variable layF'+str(nl+1)+' equal f_freeze2[2]+c_layer'+str(nl+1)+'prop[2]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')
                elif vdirec == 'x':
                    f.write('variable layF'+str(nl+1)+' equal f_freeze2[1]+c_layer'+str(nl+1)+'prop[1]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')
            elif nl == nlayer-2-1:
                if vdirec == 'y':
                    f.write('variable layF'+str(nl+1)+' equal c_layer'+str(nl)+'prop[2]+c_layer'+str(nl+1)+'prop[2]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')
                elif vdirec == 'x':
                    f.write('variable layF'+str(nl+1)+' equal c_layer'+str(nl)+'prop[1]+c_layer'+str(nl+1)+'prop[1]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')
            else:
                if vdirec == 'y':
                    f.write('variable layF'+str(nl+1)+' equal c_layer'+str(nl)+'prop[2]+c_layer'+str(nl+1)+'prop[2]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')
                elif vdirec == 'x':
                    f.write('variable layF'+str(nl+1)+' equal c_layer'+str(nl)+'prop[1]+c_layer'+str(nl+1)+'prop[1]'+'\n')
                    f.write('variable layFric'+str(nl+1)+' equal v_layF'+str(nl+1)+'/v_Pz0'+'\n')

    if nlayer != 0:
        f.write('fix thmsdupper all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_uppermsd[1] c_uppermsd[2] c_uppermsd[3] file msdupper.profile'+'\n')
        f.write('fix thmsdlower all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_lowermsd[1] c_lowermsd[2] c_lowermsd[3] file msdlower.profile'+'\n')
        f.write('fix thpeupper all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_upperpe file peupper.profile'+'\n')
        f.write('fix thpelower all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_lowerpe file pelower.profile'+'\n')
        #f.write('fix thpropupper all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_upperprop[1] c_upperprop[2] c_upperprop[3] c_upperprop[4] c_upperprop[5] c_upperprop[6] c_upperMomaxz c_upperMominz file propupper.profile'+'\n')
        #f.write('fix thproplower all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_lowerprop[1] c_lowerprop[2] c_lowerprop[3] c_lowerprop[4] c_lowerprop[5] c_lowerprop[6] c_lowerMomaxz c_lowerMominz file proplower.profile'+'\n')
        f.write('fix thpropupper all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_upperMomaxz c_upperMominz file propupper.profile'+'\n')
        f.write('fix thproplower all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_lowerMomaxz c_lowerMominz file proplower.profile'+'\n')
        for nl in range(nlayer-2):
            if nl+1 == 13:
                break
            f.write('fix thmsdlayer'+str(nl+1)+' all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_layer'+str(nl+1)+'msd[1] c_layer'+str(nl+1)+'msd[2] c_layer'+str(nl+1)+'msd[3] file msdlayer'+str(nl+1)+'.profile'+'\n')
            f.write('fix thpelayer'+str(nl+1)+' all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_layer'+str(nl+1)+'pe file pelayer'+str(nl+1)+'.profile'+'\n')
            #f.write('fix thproplayer'+str(nl+1)+' all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_layer'+str(nl+1)+'prop[1] c_layer'+str(nl+1)+'prop[2] c_layer'+str(nl+1)+'prop[3] c_layer'+str(nl+1)+'prop[4] c_layer'+str(nl+1)+'prop[5] c_layer'+str(nl+1)+'prop[6] c_layerMo'+str(nl+1)+'maxz c_layerMo'+str(nl+1)+'minz file proplayer'+str(nl+1)+'.profile'+'\n')
            f.write('fix thproplayer'+str(nl+1)+' all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' c_layerMo'+str(nl+1)+'maxz c_layerMo'+str(nl+1)+'minz file proplayer'+str(nl+1)+'.profile'+'\n')
            f.write('fix thlayfric'+str(nl+1)+' all ave/time 1 '+str(int(md_step/nave))+' '+str(int(md_step/nave))+' v_layFric'+str(nl+1)+' file friclay'+str(nl+1)+'.profile'+'\n')
            #if nl == 0:
            #    f.write('fix thuptol'+str(nl+1)+' all ave/time '+str(dmp_step)+' '+str(dmp_step)+' '+str(md_step)+' c_uptol'+str(nl+1)+' c_uptol'+str(nl+1)+'[1] c_uptol'+str(nl+1)+'[2] c_uptol'+str(nl+1)+'[3] file uptol'+str(nl+1)+'.profile'+'\n')
            #elif nl == nlayer-2-1:
            #    f.write('fix thl'+str(nl+1)+'tolow all ave/time '+str(dmp_step)+' '+str(dmp_step)+' '+str(md_step)+' c_l'+str(nl+1)+'tolow c_l'+str(nl+1)+'tolow[1] c_l'+str(nl+1)+'tolow[2] c_l'+str(nl+1)+'tolow[3] file l'+str(nl+1)+'tolow.profile'+'\n')
            #else:
            #    f.write('fix thl'+str(nl)+'tol'+str(nl+1)+' all ave/time '+str(dmp_step)+' '+str(dmp_step)+' '+str(md_step)+' c_l'+str(nl)+'tol'+str(nl+1)+' c_l'+str(nl)+'tol'+str(nl+1)+'[1] c_l'+str(nl)+'tol'+str(nl+1)+'[2] c_l'+str(nl)+'tol'+str(nl+1)+'[3] file l'+str(nl)+'tol'+str(nl+1)+'.profile'+'\n')
    if nlayer != 0:
        f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma f_freeze2[1] f_freeze2[2] f_freeze2[3] v_Pz0 v_Fric')
        for nl in range(nlayer-2):
            if nl+1 == 13:
                break
            f.write(' v_layFric'+str(nl+1))
        f.write('\n')
    else:
        f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma f_freeze2[1] f_freeze2[2] f_freeze2[3] v_Pz0 v_Fric'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('fix fxnvt middle nvt temp '+str(T)+' '+str(T)+' $(dt*100.0)'+'\n')
    f.write('run '+str(md_step)+'\n')
    f.close()

def mk_npt_melt_input_deepmd(cell,dt,dmp_step,thermo_step,eq_step,melt_step,md_step,Teq,Tmelt,T,P,z0,seed):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    mk_lammpsdata(cell,False,force_skew=True)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('region freeze block INF INF INF INF INF '+str(z0)+'\n')
    f.write('group lower region freeze'+'\n')
    f.write('group upper subtract all lower'+'\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(Teq)+' '+str(seed)+'\n')
    f.write('fix fmom all momentum 1 linear 1 1 1'+'\n')
    f.write('fix fxnpt all npt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('run '+str(eq_step)+'\n')
    f.write('unfix fxnpt'+'\n')
    f.write('velocity upper scale '+str(Tmelt)+'\n')
    f.write('fix fxnvt upper nvt temp '+str(Tmelt)+' '+str(Tmelt)+' $(dt*100.0)'+'\n')
    f.write('run '+str(melt_step)+'\n')
    f.write('unfix fxnvt'+'\n')
    f.write('write_data result0.data'+'\n')
    f.write('velocity all scale '+str(T)+'\n')
    f.write('fix th all ave/time 1 '+str(int(md_step/2))+' '+str(md_step)+' v_mytemp v_myenthalpy v_mydensity v_myvol v_mycella v_mycellb v_mycellc v_mycellalpha v_mycellbeta v_mycellgamma file th'+str(T)+'.profile'+'\n')
    f.write('fix fxnpt all npt temp '+str(T)+' '+str(T)+' $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('run '+str(md_step)+'\n')
    f.write('unfix fxnpt'+'\n')
    f.write('unfix th'+'\n')
    f.close()

def mk_npt_input_deepmd(cell,dt,dmp_step,thermo_step,eq_step,Teq,P,seed,mol=False):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    mk_lammpsdata(cell,False,force_skew=True,mol=mol)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(Teq)+' '+str(seed)+'\n')
    f.write('fix fmom all momentum 1 linear 1 1 1'+'\n')
    f.write('fix fxnpt all npt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('run '+str(eq_step)+'\n')
    f.write('write_data result.data'+'\n')
    f.close()

def mk_npt_input_dpdata(cell,dt,thermo_step,eq_step,md_step,nout,Teq,P,seed,mol=False):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    mk_lammpsdata(cell,False,force_skew=True,mol=mol)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('variable mype equal pe'+'\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    #f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id type element x y z ix iy iz'+'\n')
    #f.write('dump_modify dmp element')
    #for symbol in symbols:
    #    f.write(' '+str(symbol))
    #f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(Teq)+' '+str(seed)+'\n')
    #f.write('velocity all create 300 '+str(seed)+'\n')
    f.write('fix fmom all momentum 1 linear 1 1 1'+'\n')
    f.write('fix fxnpt all npt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    #f.write('fix fxnpt all npt temp 300 300 $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('run '+str(eq_step)+'\n')
    f.write('unfix fxnpt'+'\n')
    f.write('label loop'+'\n')
    f.write('  variable a loop '+str(nout)+'\n')
    f.write('  reset_timestep 0'+'\n')
    f.write('  dump myDump all custom '+str(md_step)+' dump${a}.lammpstrj.* id type fx fy fz'+'\n')
    f.write('  fix fxnpt all npt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('  run '+str(md_step)+'\n')
    f.write('  print "step ${a} Potential Energy: ${mype}"'+'\n')
    f.write('  write_data result${a}.data'+'\n')
    f.write('  undump myDump'+'\n')
    f.write('  unfix fxnpt'+'\n')
    f.write('  next a'+'\n')
    f.write('jump SELF loop'+'\n')
    f.close()

def mk_nvt_input_dpdata(cell,dt,thermo_step,eq_step,md_step,nout,Teq,seed,mol=False):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    mk_lammpsdata(cell,False,force_skew=True,mol=mol)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    f.write('boundary p p p'+'\n')
    f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('variable mype equal pe'+'\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    #f.write('dump dmp all custom '+str(dmp_step)+' traj_npt.lammpstrj id type element x y z ix iy iz'+'\n')
    #f.write('dump_modify dmp element')
    #for symbol in symbols:
    #    f.write(' '+str(symbol))
    #f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('timestep '+str(dt)+'\n')
    f.write('velocity all create '+str(Teq)+' '+str(seed)+'\n')
    #f.write('velocity all create 300 '+str(seed)+'\n')
    f.write('fix fmom all momentum 1 linear 1 1 1'+'\n')
    f.write('fix fxnvt all nvt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0)'+'\n')
    #f.write('fix fxnpt all npt temp 300 300 $(dt*100.0) '+'tri '+str(P)+' '+str(P)+' $(dt*1000.0)'+'\n')
    f.write('run '+str(eq_step)+'\n')
    f.write('unfix fxnvt'+'\n')
    f.write('label loop'+'\n')
    f.write('  variable a loop '+str(nout)+'\n')
    f.write('  reset_timestep 0'+'\n')
    f.write('  dump myDump all custom '+str(md_step)+' dump${a}.lammpstrj.* id type fx fy fz'+'\n')
    f.write('  fix fxnvt all nvt temp '+str(Teq)+' '+str(Teq)+' $(dt*100.0)'+'\n')
    f.write('  run '+str(md_step)+'\n')
    f.write('  print "step ${a} Potential Energy: ${mype}"'+'\n')
    f.write('  write_data result${a}.data'+'\n')
    f.write('  undump myDump'+'\n')
    f.write('  unfix fxnvt'+'\n')
    f.write('  next a'+'\n')
    f.write('jump SELF loop'+'\n')
    f.close()

def mk_opt_input_deepmd(cell,dmp_step,thermo_step,mol=False,fixlay=None):
    from AtomicVirtuaLab.io import cell2atomlist, mk_lammpsdata
    from ase import Atom
    mk_lammpsdata(cell,False,force_skew=True,mol=mol)
    symbols = cell2atomlist(cell)
    f = open('lammps.lmp','w')
    f.write('units metal'+'\n')
    if fixlay == None:
        f.write('boundary p p p'+'\n')
    else:
        f.write('boundary p p f'+'\n')
    if mol:
        f.write('atom_style full'+'\n')
    else:
        f.write('atom_style atomic'+'\n')
    f.write('read_data lammps.data'+'\n')
    i = 1
    for symbol in symbols:
        a = Atom(symbol)
        f.write('mass '+str(i)+' '+str(a.mass)+'\n')
        i = i + 1
    if fixlay != None:
        f.write('region lower block INF INF INF INF INF '+str(fixlay)+'\n')
        f.write('group lower region lower'+'\n')
        f.write('fix freeze lower setforce 0.0 0.0 0.0'+'\n')
    f.write('pair_style deepmd graph.pb'+'\n')
    f.write('pair_coeff * * ')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('variable mytemp equal temp'+'\n')
    f.write('variable myenthalpy equal enthalpy'+'\n')
    f.write('variable mydensity equal density'+'\n')
    f.write('variable myvol equal vol'+'\n')
    f.write('variable mycella equal cella'+'\n')
    f.write('variable mycellb equal cellb'+'\n')
    f.write('variable mycellc equal cellc'+'\n')
    f.write('variable mycellalpha equal cellalpha'+'\n')
    f.write('variable mycellbeta equal cellbeta'+'\n')
    f.write('variable mycellgamma equal cellgamma'+'\n')
    f.write('dump dmp all custom '+str(dmp_step)+' traj_min.lammpstrj id type element x y z ix iy iz'+'\n')
    f.write('dump_modify dmp element')
    for symbol in symbols:
        f.write(' '+str(symbol))
    f.write('\n')
    f.write('thermo_style custom step etotal enthalpy pe ke temp press vol density cella cellb cellc cellalpha cellbeta cellgamma'+'\n')
    f.write('thermo '+str(thermo_step)+'\n')
    f.write('\n')
    f.write('minimize 0.0 1.0e-12 10000 100000'+'\n')
    f.write('write_data result.data'+'\n')
    f.write('variable mype equal pe'+'\n')
    f.write('run 0'+'\n')
    f.write('print "Potential Energy: ${mype}"'+'\n')
    f.close()

def mk_qeqfile(symbols,unit='metal'):
    from pymatgen.core.periodic_table import Element
    f = open('my_qeq','w')
    if unit == 'metal':
        f.write('# UNITS: metal'+'\n')
    elif unit == 'real':
        f.write('# UNITS: real'+'\n')
    i = 1
    for symbol in symbols:
        IP = Element(symbol).data['Ionization energies'][0]
        EA = Element(symbol).data['Electron affinity']
        CHI = (IP+EA)/2.0
        ETA = (IP-EA)/2.0
        if unit == 'real':
            CHI = CHI*23.06000
            ETA = ETA*23.06000
        f.write(str(i)+' '+str(CHI)+' '+str(ETA)+' '+'0.0 0.0 0.0'+'\n')
        i = i + 1
    f.close()

def set_uff_lj(atom):
    scale = 1.0/(2.0)**(1.0/6.0)
    kcal2eV = 4.336454e-2
    ljparams={
        'H':{
            'sigma'  : 2.886*scale,
            'epsilon': 0.044*kcal2eV
        },
        'He':{
            'sigma'  : 2.362*scale,
            'epsilon': 0.056*kcal2eV
        },
        'Li':{
            'sigma'  : 2.451*scale,
            'epsilon': 0.025*kcal2eV
        },
        'Be':{
            'sigma'  : 2.745*scale,
            'epsilon': 0.085*kcal2eV
        },
        'B':{
            'sigma'  : 4.083*scale,
            'epsilon': 0.180*kcal2eV
        },
        'C':{
            'sigma'  : 3.851*scale,
            'epsilon': 0.105*kcal2eV
        },
        'N':{
            'sigma'  : 3.660*scale,
            'epsilon': 0.069*kcal2eV
        },
        'O':{
            'sigma'  : 3.500*scale,
            'epsilon': 0.060*kcal2eV
        },
        'F':{
            'sigma'  : 3.364*scale,
            'epsilon': 0.050*kcal2eV
        },
        'Ne':{
            'sigma'  : 3.243*scale,
            'epsilon': 0.042*kcal2eV
        },
        'Na':{
            'sigma'  : 2.983*scale,
            'epsilon': 0.030*kcal2eV
        },
        'Mg':{
            'sigma'  : 3.021*scale,
            'epsilon': 0.111*kcal2eV
        },
        'Al':{
            'sigma'  : 4.499*scale,
            'epsilon': 0.505*kcal2eV
        },
        'Si':{
            'sigma'  : 4.295*scale,
            'epsilon': 0.402*kcal2eV
        },
        'P':{
            'sigma'  : 4.147*scale,
            'epsilon': 0.305*kcal2eV
        },
        'S':{
            'sigma'  : 4.035*scale,
            'epsilon': 0.274*kcal2eV
        },
        'Cl':{
            'sigma'  : 3.947*scale,
            'epsilon': 0.227*kcal2eV
        },
        'Ar':{
            'sigma'  : 3.868*scale,
            'epsilon': 0.185*kcal2eV
        },
        'K':{
            'sigma'  : 3.812*scale,
            'epsilon': 0.035*kcal2eV
        },
        'Ca':{
            'sigma'  : 3.399*scale,
            'epsilon': 0.238*kcal2eV
        },
        'Sc':{
            'sigma'  : 3.295*scale,
            'epsilon': 0.019*kcal2eV
        },
        'Ti':{
            'sigma'  : 3.175*scale,
            'epsilon': 0.017*kcal2eV
        },
        'V':{
            'sigma'  : 3.144*scale,
            'epsilon': 0.016*kcal2eV
        },
        'Cr':{
            'sigma'  : 3.023*scale,
            'epsilon': 0.015*kcal2eV
        },
        'Mn':{
            'sigma'  : 2.961*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Fe':{
            'sigma'  : 2.912*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Co':{
            'sigma'  : 2.872*scale,
            'epsilon': 0.014*kcal2eV
        },
        'Ni':{
            'sigma'  : 2.834*scale,
            'epsilon': 0.015*kcal2eV
        },
        'Cu':{
            'sigma'  : 3.495*scale,
            'epsilon': 0.005*kcal2eV
        },
        'Zn':{
            'sigma'  : 2.763*scale,
            'epsilon': 0.124*kcal2eV
        },
        'Ga':{
            'sigma'  : 4.383*scale,
            'epsilon': 0.415*kcal2eV
        },
        'Ge':{
            'sigma'  : 4.280*scale,
            'epsilon': 0.379*kcal2eV
        },
        'As':{
            'sigma'  : 4.230*scale,
            'epsilon': 0.309*kcal2eV
        },
        'Se':{
            'sigma'  : 4.205*scale,
            'epsilon': 0.291*kcal2eV
        },
        'Br':{
            'sigma'  : 4.189*scale,
            'epsilon': 0.251*kcal2eV
        },
        'Kr':{
            'sigma'  : 4.141*scale,
            'epsilon': 0.220*kcal2eV
        },
        'Rb':{
            'sigma'  : 4.114*scale,
            'epsilon': 0.040*kcal2eV
        },
        'Sr':{
            'sigma'  : 3.641*scale,
            'epsilon': 0.235*kcal2eV
        },
        'Y':{
            'sigma'  : 3.345*scale,
            'epsilon': 0.072*kcal2eV
        },
        'Zr':{
            'sigma'  : 3.124*scale,
            'epsilon': 0.069*kcal2eV
        },
        'Nb':{
            'sigma'  : 3.165*scale,
            'epsilon': 0.059*kcal2eV
        },
        'Mo':{
            'sigma'  : 3.052*scale,
            'epsilon': 0.056*kcal2eV
        },
        'Tc':{
            'sigma'  : 2.998*scale,
            'epsilon': 0.048*kcal2eV
        },
        'Ru':{
            'sigma'  : 2.963*scale,
            'epsilon': 0.056*kcal2eV
        },
        'Rh':{
            'sigma'  : 2.929*scale,
            'epsilon': 0.053*kcal2eV
        },
        'Pd':{
            'sigma'  : 2.899*scale,
            'epsilon': 0.048*kcal2eV
        },
        'Ag':{
            'sigma'  : 3.148*scale,
            'epsilon': 0.036*kcal2eV
        },
        'Cd':{
            'sigma'  : 2.848*scale,
            'epsilon': 0.228*kcal2eV
        },
        'In':{
            'sigma'  : 4.463*scale,
            'epsilon': 0.599*kcal2eV
        },
        'Sn':{
            'sigma'  : 4.392*scale,
            'epsilon': 0.567*kcal2eV
        },
        'Sb':{
            'sigma'  : 4.420*scale,
            'epsilon': 0.449*kcal2eV
        },
        'Te':{
            'sigma'  : 4.470*scale,
            'epsilon': 0.398*kcal2eV
        },
        'I':{
            'sigma'  : 4.500*scale,
            'epsilon': 0.339*kcal2eV
        },
        'Xe':{
            'sigma'  : 4.404*scale,
            'epsilon': 0.332*kcal2eV
        },
        'Cs':{
            'sigma'  : 4.517*scale,
            'epsilon': 0.045*kcal2eV
        },
        'Ba':{
            'sigma'  : 3.703*scale,
            'epsilon': 0.364*kcal2eV
        },
        'La':{
            'sigma'  : 3.522*scale,
            'epsilon': 0.017*kcal2eV
        },
        'Ce':{
            'sigma'  : 3.556*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Pr':{
            'sigma'  : 3.606*scale,
            'epsilon': 0.010*kcal2eV
        },
        'Nd':{
            'sigma'  : 3.575*scale,
            'epsilon': 0.010*kcal2eV
        },
        'Pm':{
            'sigma'  : 3.547*scale,
            'epsilon': 0.009*kcal2eV
        },
        'Sm':{
            'sigma'  : 3.520*scale,
            'epsilon': 0.008*kcal2eV
        },
        'Eu':{
            'sigma'  : 3.493*scale,
            'epsilon': 0.008*kcal2eV
        },
        'Gd':{
            'sigma'  : 3.368*scale,
            'epsilon': 0.009*kcal2eV
        },
        'Tb':{
            'sigma'  : 3.451*scale,
            'epsilon': 0.007*kcal2eV
        },
        'Dy':{
            'sigma'  : 3.428*scale,
            'epsilon': 0.007*kcal2eV
        },
        'Ho':{
            'sigma'  : 3.409*scale,
            'epsilon': 0.007*kcal2eV
        },
        'Er':{
            'sigma'  : 3.391*scale,
            'epsilon': 0.007*kcal2eV
        },
        'Tm':{
            'sigma'  : 3.374*scale,
            'epsilon': 0.006*kcal2eV
        },
        'Yb':{
            'sigma'  : 3.355*scale,
            'epsilon': 0.228*kcal2eV
        },
        'Lu':{
            'sigma'  : 3.640*scale,
            'epsilon': 0.041*kcal2eV
        },
        'Hf':{
            'sigma'  : 3.141*scale,
            'epsilon': 0.072*kcal2eV
        },
        'Ta':{
            'sigma'  : 3.170*scale,
            'epsilon': 0.081*kcal2eV
        },
        'W':{
            'sigma'  : 3.069*scale,
            'epsilon': 0.067*kcal2eV
        },
        'Re':{
            'sigma'  : 2.954*scale,
            'epsilon': 0.066*kcal2eV
        },
        'Os':{
            'sigma'  : 3.120*scale,
            'epsilon': 0.037*kcal2eV
        },
        'Ir':{
            'sigma'  : 2.840*scale,
            'epsilon': 0.073*kcal2eV
        },
        'Pt':{
            'sigma'  : 2.754*scale,
            'epsilon': 0.080*kcal2eV
        },
        'Au':{
            'sigma'  : 3.293*scale,
            'epsilon': 0.039*kcal2eV
        },
        'Hg':{
            'sigma'  : 2.705*scale,
            'epsilon': 0.385*kcal2eV
        },
        'Tl':{
            'sigma'  : 4.347*scale,
            'epsilon': 0.680*kcal2eV
        },
        'Pb':{
            'sigma'  : 4.297*scale,
            'epsilon': 0.663*kcal2eV
        },
        'Bi':{
            'sigma'  : 4.370*scale,
            'epsilon': 0.518*kcal2eV
        },
        'Po':{
            'sigma'  : 4.709*scale,
            'epsilon': 0.325*kcal2eV
        },
        'At':{
            'sigma'  : 4.750*scale,
            'epsilon': 0.284*kcal2eV
        },
        'Rn':{
            'sigma'  : 4.765*scale,
            'epsilon': 0.248*kcal2eV
        },
        'Fr':{
            'sigma'  : 4.900*scale,
            'epsilon': 0.050*kcal2eV
        },
        'Ra':{
            'sigma'  : 3.677*scale,
            'epsilon': 0.404*kcal2eV
        },
        'Ac':{
            'sigma'  : 3.478*scale,
            'epsilon': 0.033*kcal2eV
        },
        'Th':{
            'sigma'  : 3.396*scale,
            'epsilon': 0.026*kcal2eV
        },
        'Pa':{
            'sigma'  : 3.424*scale,
            'epsilon': 0.022*kcal2eV
        },
        'U':{
            'sigma'  : 3.395*scale,
            'epsilon': 0.022*kcal2eV
        },
        'Np':{
            'sigma'  : 3.424*scale,
            'epsilon': 0.019*kcal2eV
        },
        'Pu':{
            'sigma'  : 3.424*scale,
            'epsilon': 0.016*kcal2eV
        },
        'Am':{
            'sigma'  : 3.381*scale,
            'epsilon': 0.014*kcal2eV
        },
        'Cm':{
            'sigma'  : 3.326*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Bk':{
            'sigma'  : 3.339*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Cf':{
            'sigma'  : 3.313*scale,
            'epsilon': 0.013*kcal2eV
        },
        'Es':{
            'sigma'  : 3.299*scale,
            'epsilon': 0.012*kcal2eV
        },
        'Fm':{
            'sigma'  : 3.286*scale,
            'epsilon': 0.012*kcal2eV
        },
        'Md':{
            'sigma'  : 3.274*scale,
            'epsilon': 0.011*kcal2eV
        },
        'No':{
            'sigma'  : 3.248*scale,
            'epsilon': 0.011*kcal2eV
        },
        'Lw':{
            'sigma'  : 3.236*scale,
            'epsilon': 0.011*kcal2eV
        }
    }
    return ljparams[atom]['sigma'], ljparams[atom]['epsilon']
