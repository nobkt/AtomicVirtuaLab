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
    ax.get_legend().remove()
    #ax.legend(loc="upper right")
    plt.xlabel("Wave vector", fontsize=16)
    plt.ylabel("E-Efermi(eV)", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.savefig('vasp_band.png')
    #plt.show()
    plt.clf()

def plot_vasp_dos(fvaspxml):
    from pymatgen.io.vasp import Vasprun
    from pymatgen.electronic_structure.plotter import DosPlotter
    from matplotlib import pyplot as plt
    import sys
    v = Vasprun(fvaspxml)
    tdos = v.tdos
    x = list(tdos.energies-tdos.efermi)
    for spin in tdos.densities:
        if spin == spin.up:
            y = list(tdos.densities[spin])
            plt.plot(x,y,label='Total DOS(UP)', linewidth = 0.8)
        elif spin == spin.down:
            y = list(tdos.densities[spin])
            plt.plot(x,y,label='Total DOS(DOWN)', linewidth = 0.8)
    plt.xlabel('E-Efermi(eV)', fontsize=16)
    plt.ylabel('DOS', fontsize=16)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black", linewidth = 0.6)
    plt.xlim(-5.0,5.0)
    plt.ylim(0,20.0)
    plt.tick_params(labelsize=14)
    plt.legend(loc='upper right',fontsize=9)
    plt.savefig('vasp_total_dos.png',dpi=500)
    #plt.show()
    plt.clf()
    
    cdos = v.complete_dos
    element_dos = cdos.get_element_dos()
    for el in element_dos:
        el_spd = cdos.get_element_spd_dos(el)
        for spd in el_spd:
            x = list(el_spd[spd].energies-el_spd[spd].efermi)
            for spin in el_spd[spd].densities:
                if spin == spin.up:
                    y = list(el_spd[spd].densities[spin])
                    plt.plot(x,y,label=str(el)+'_'+str(spd)+'(UP)', linewidth = 0.8)
                elif spin == spin.down:
                    y = list(el_spd[spd].densities[spin])
                    plt.plot(x,y,label=str(el)+'_'+str(spd)+'(DOWN)', linewidth = 0.8)
    plt.xlabel('E-Efermi(eV)', fontsize=16)
    plt.ylabel('DOS', fontsize=16)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black", linewidth = 0.6)
    plt.xlim(-5.0,5.0)
    plt.ylim(0,20.0)
    plt.tick_params(labelsize=14)
    plt.legend(loc='upper right',fontsize=9)
    plt.savefig('vasp_pdos.png',dpi=500)
    #plt.show()
    plt.clf()
    
    pdos = cdos.pdos
    merged = {}
    for site in pdos:
        orbitals = pdos[site]
        site0 = str(site)
        atom = site0[site0.index(']') + 2:]
        for orbital in orbitals:
            orb_dos = cdos.get_site_orbital_dos(site,orbital)
            for spin in orb_dos.densities:
                if spin == spin.up:
                    mkey = str(atom)+'_'+str(orbital)+'(UP)'
                elif spin == spin.down:
                    mkey = str(atom)+'_'+str(orbital)+'(DOWN)'
                if mkey not in merged.keys():
                    merged[mkey] = orb_dos
                else:
                    merged[mkey] += orb_dos
    for mkey in merged:
        orb_dos = merged[mkey]
        x = list(orb_dos.energies-orb_dos.efermi)
        for spin in orb_dos.densities:
            y = list(orb_dos.densities[spin])
            plt.plot(x,y,label=mkey, linewidth = 0.8)
        plt.xlabel('E-Efermi(eV)', fontsize=16)
    plt.ylabel('DOS', fontsize=16)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black", linewidth = 0.6)
    plt.xlim(-5.0,5.0)
    plt.ylim(0,20.0)
    plt.tick_params(labelsize=14)
    plt.legend(loc='upper right',fontsize=9)
    plt.savefig('vasp_orbital_dos.png',dpi=500)
    #plt.show()
    plt.clf()

    """
    cdos = v.complete_dos
    pdos = cdos.pdos
    plotter = DosPlotter()
    plotter2 = DosPlotter()
    element_dos = cdos.get_element_dos()
    #for el in element_dos:
    #    el_spd = cdos.get_element_spd_dos(el)
    #    for spd in el_spd:
    #        print(el,spd,el_spd[spd])
    merged = {}
    for site in pdos:
        orbitals = pdos[site]
        site0 = str(site)
        atom = site0[site0.index(']') + 2:]
        for orbital in orbitals:
            #print(site)
            #print(orbital)
            #print(atom)
            orb_dos = cdos.get_site_orbital_dos(site,orbital)
            mkey = str(atom)+'_'+str(orbital)
            if mkey not in merged.keys():
                merged[mkey] = orb_dos
            else:
                merged[mkey] += orb_dos
    for mkey in merged:
        orb_dos = merged[mkey]
        plotter2.add_dos(mkey,orb_dos)
    plotter.add_dos_dict(element_dos)
    plotter.add_dos("Total DOS", tdos)
    plt = plotter2.get_plot(xlim=[-5.0,5.0],ylim=[0,20.0])
    plt.xlabel("xlabel", fontsize=16)
    plt.ylabel("ylabel", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.legend(loc='upper right')
    plt.savefig('dos.png',dpi=500)
    plt.show()
    """




