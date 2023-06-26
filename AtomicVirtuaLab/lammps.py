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
    
def mk_npt_input_fr_moltemplate(symbols,minimize,dt,dmp_step,thermo_step,md_step,T,P,dof_type,seed,restart):
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

def mk_qeqfile(symbols):
    from pymatgen.core.periodic_table import Element
    f = open('my_qeq','w')
    f.write('# UNITS: metal'+'\n')
    i = 1
    for symbol in symbols:
        IP = Element(symbol).data['Ionization energies'][0]
        EA = Element(symbol).data['Electron affinity']
        CHI = (IP+EA)/2.0
        ETA = (IP-EA)/2.0
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
