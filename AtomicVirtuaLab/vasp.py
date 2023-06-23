def mk_poscar(cell):
    from ase.io import write
    cell.write('POSCAR',format='vasp',sort=True)

def mk_potcar(cell,potdir):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    import os
    symbols = cell2atomlist(cell)
    command = 'cat '
    for symbol in symbols:
        command = command + potdir+'/'+symbol+'/POTCAR '
    command = command + '> POTCAR'
    os.system(command)

def mk_npt_incar(cell,T,P,nsw,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ISMEAR = 0'+'\n')
    f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('IBRION = 0'+'\n')
    f.write('MDALGO = 3'+'\n')
    f.write('ISIF = 3'+'\n')
    f.write('TEBEG = '+str(T)+'\n')
    f.write('TEEND = '+str(T)+'\n')
    f.write('NSW = '+str(nsw)+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('POTIM = 4.0'+'\n')
    f.write('LANGEVIN_GAMMA =')
    for symbol in symbols:
        f.write(' 10.0')
    f.write('\n')
    f.write('LANGEVIN_GAMMA_L = 1.0'+'\n')
    f.write('PSTRESS = '+str(P)+'\n')
    if fscan:
        f.write('GGA = '+str(GGA)+'\n')
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_nvt_incar(cell,T,nsw,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ISMEAR = 0'+'\n')
    f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('IBRION = 0'+'\n')
    f.write('MDALGO = 3'+'\n')
    f.write('ISIF = 2'+'\n')
    f.write('TEBEG = '+str(T)+'\n')
    f.write('TEEND = '+str(T)+'\n')
    f.write('NSW = '+str(nsw)+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('POTIM = 2.0'+'\n')
    f.write('LANGEVIN_GAMMA =')
    for symbol in symbols:
        f.write(' 10.0')
    f.write('\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_opt_incar(cell,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ISMEAR = 0'+'\n')
    f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    f.write('ISIF = 2'+'\n')
    f.write('IBRION = 2'+'\n')
    f.write('NSW = 10000'+'\n')
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_cellopt_incar(cell,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ISMEAR = 0'+'\n')
    f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    f.write('ISIF = 3'+'\n')
    f.write('IBRION = 2'+'\n')
    f.write('NSW = 10000'+'\n')
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_scf_incar(cell,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ISMEAR = 0'+'\n')
    f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_dos_incar(cell,GGA='PE',fscan=False,options=[],ispin=1):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ICHARG = 11'+'\n')
    f.write('ISMEAR = -5'+'\n')
    f.write('LORBIT = 11'+'\n')
    #f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_band_incar(cell,GGA='PE',fscan=False,options=[],ispin=1,PREC='Normal'):
    from AtomicVirtuaLab.io import cell2atomlist
    import re
    symbols = cell2atomlist(cell)
    f = open('INCAR','w')
    f.write('NCORE = 4'+'\n')
    f.write('ICHARG = 11'+'\n')
    f.write('ISMEAR = 0'+'\n')
    #f.write('LREAL = Auto'+'\n')
    f.write('ALGO = Normal'+'\n')
    f.write('PREC = Normal'+'\n')
    f.write('LASPH = .TRUE.'+'\n')
    f.write('ISYM = 0'+'\n')
    f.write('NELM = 10000'+'\n')
    f.write('ISPIN = '+str(ispin)+"\n")
    if len(options) != 0:
        for option in options:
                f.write(option+'\n')
    if fscan:
        f.write('METAGGA=SCAN'+'\n')
        f.write('AMIX=0.1'+'\n')
        f.write('GGA = '+str(GGA)+'\n')
    else:
        f.write('GGA = '+str(GGA)+'\n')
    f.close()

def mk_kpoints(kx,ky,kz,band=False):
    if band == False:
        f = open('KPOINTS','w')
        f.write('K-Points'+'\n')
        f.write(' 0'+'\n')
        f.write('Gamma'+'\n')
        f.write(' '+str(kx)+'  '+str(ky)+'  '+str(kz)+'\n')
        f.write(' 0  0  0'+'\n')
        f.close()
    elif band == True:
        from pymatgen.io.vasp.inputs import Kpoints
        from pymatgen.core import Structure
        from pymatgen.symmetry.bandstructure import HighSymmKpath
        struct = Structure.from_file("POSCAR")
        kpath = HighSymmKpath(struct)
        kpts = Kpoints.automatic_linemode(divisions=40,ibz=kpath)
        kpts.write_file("KPOINTS")

def plot_vasp_band(fvaspxml):
    from pymatgen.io.vasp.outputs import Vasprun
    from pymatgen.electronic_structure.plotter import BSPlotter
    vaspout = Vasprun(fvaspxml)
    bandstr = vaspout.get_band_structure(line_mode=True)
    plt = BSPlotter(bandstr).get_plot(ylim=[-10,10])
    ax = plt.gca()
    xlim = ax.get_xlim()
    ax.hlines(0,xlim[0],xlim[1], linestyles="dashed", color="black", linewidth = 0.6)
    ax.legend(loc="upper left")
    plt.xlabel("xlabel", fontsize=16)
    plt.ylabel("ylabel", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.savefig('band.png')
    plt.show()

def plot_vasp_dos(fvaspxml):
    from pymatgen.io.vasp import Vasprun
    from pymatgen.electronic_structure.plotter import DosPlotter
    v = Vasprun(fvaspxml)
    tdos = v.tdos
    cdos = v.complete_dos
    pdos = cdos.pdos
    plotter = DosPlotter()
    element_dos = cdos.get_element_dos()
    plotter.add_dos("Total DOS", tdos)
    plotter.add_dos_dict(element_dos)
    merged = {}
    for site in pdos:
        orbitals = pdos[site]
        site0 = str(site)
        atom = site0[site0.index(']') + 2:]
        for orbital in orbitals:
            print(site)
            print(orbital)
            print(atom)
            orb_dos = cdos.get_site_orbital_dos(site,orbital)
            mkey = str(atom)+'_'+str(orbital)
            if mkey not in merged.keys():
                merged[mkey] = orb_dos
            else:
                merged[mkey] += orb_dos
    for mkey in merged:
        orb_dos = merged[mkey]
        plotter.add_dos(mkey,orb_dos)
    plt = plotter.get_plot()
    plt.xlabel("xlabel", fontsize=16)
    plt.ylabel("ylabel", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.savefig('band.png')
    plt.legend(loc='upper left')
    plt.xlim(-5.0,5.0)
    plt.ylim(0,20.0)
    plt.show()



