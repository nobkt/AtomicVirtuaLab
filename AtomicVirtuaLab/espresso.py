def mk_qe_input_scf(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'scf',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : True,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_scf.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_dos(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'verbosity'        : 'high',\
        'calculation'      : 'nscf',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        #'tstress'          : True,\
        #'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'occupations'      : 'tetrahedra',\
        #'occupations'      : 'smearing',\
        #'smearing'         : 'gaussian',\
        #'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_nscf_dos.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)
    f = open('qe_dos.pwi','w')
    f.write('&DOS'+'\n')
    f.write('outdir = "./outdir"'+'\n')
    #f.write('degauss =  0.01'+'\n')
    f.write('deltae  =  0.01'+'\n')
    #f.write('emax    =  50.0'+'\n')
    #f.write('emin    = -50.0'+'\n')
    #f.write('ngauss  =  0'+'\n')
    f.write('/'+'\n')
    f.close()
    f = open('qe_projwfc.pwi','w')
    f.write('&PROJWFC'+'\n')
    f.write('outdir = "./outdir"'+'\n')
    #f.write('degauss =  0.01'+'\n')
    f.write('deltae  =  0.01'+'\n')
    #f.write('emax    =  50.0'+'\n')
    #f.write('emin    = -50.0'+'\n')
    #f.write('ngauss  =  0'+'\n')
    f.write('/'+'\n')
    f.close()

def mk_qe_input_band(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import sys
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'verbosity'        : 'high',\
        'calculation'      : 'bands',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        #'tstress'          : True,\
        #'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    lat = cell.cell.get_bravais_lattice()
    paths = lat.special_path
    npoints = len(paths) - paths.count(',')
    bandpath = {'path':paths, 'npoints':npoints-1}
    cell.write('qe_nscf_band.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=bandpath,crystal_coordinates=False)
    f = open('qe_nscf_band.pwi','r')
    lines = f.readlines()
    f.close()
    f = open('qe_nscf_band.pwi','w')
    iline = 0
    ipath = 0
    wflg1 = 0
    wflg2 = 0
    for line in lines:
        if 'K_POINTS' in line:
            wflg1 = 1
        if wflg1 == 0:
            f.write(line)
        elif wflg1 == 1 and wflg2 < 2:
            f.write(line)
            wflg2 = wflg2 + 1
        elif wflg1 == 1 and wflg2 == 2:
            if paths[ipath] == ',':
                ipath = ipath + 1
            line0 = line.split()
            npoint = 20
            if ipath+1 > len(paths)-1:
                npoint = 0
                wflg1 = 0
            elif paths[ipath+1] == ',':
                npoint = 0
            f.write(line0[0]+' '+line0[1]+' '+line0[2]+' '+str(npoint)+' !'+paths[ipath]+'\n')
            ipath = ipath + 1
        iline = iline + 1
    f.close()
    if nspin == False:
        f = open('qe_band.pwi','w')
        f.write('&BANDS'+'\n')
        f.write('outdir = "./outdir"'+'\n')
        f.write('lsym =  .FALSE.'+'\n')
        f.write('spin_component =  1'+'\n')
        f.write('/'+'\n')
        f.close()
    elif nspin == True:
        f = open('qe_band_up.pwi','w')
        f.write('&BANDS'+'\n')
        f.write('outdir = "./outdir"'+'\n')
        f.write('filband = "bands_up.out"'+'\n')
        f.write('lsym =  .FALSE.'+'\n')
        f.write('spin_component =  1'+'\n')
        f.write('/'+'\n')
        f.close()
        f = open('qe_band_down.pwi','w')
        f.write('&BANDS'+'\n')
        f.write('outdir = "./outdir"'+'\n')
        f.write('filband = "bands_down.out"'+'\n')
        f.write('lsym =  .FALSE.'+'\n')
        f.write('spin_component =  2'+'\n')
        f.write('/'+'\n')
        f.close()

def mk_qe_input_relax(cell,xc,pot,level='low',estep=1000,nstep=1000,nosym=False,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'relax',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : True,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nstep'            : nstep,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'nosym'            : nosym,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_relax.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_vcrelax(cell,xc,pot,level='low',estep=1000,nstep=1000,nosym=False,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'vc-relax',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : True,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nstep'            : nstep,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'nosym'            : nosym,\
        'cell_dofree'      : 'all',\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_vc-relax.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_nvt(cell,xc,pot,tempw,tolp,dt=0.5,level='low',estep=1000,nstep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'md',\
        'dt'               : dt*1.0e-15/4.8378e-17,\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : True,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nstep'            : nstep,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'nosym'            : True,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_,\
        'ion_dynamics'     : 'verlet',\
        'ion_temperature'  : 'rescaling',\
        'tempw'            : tempw,\
        'tolp'             : tolp,\
        'pot_extrapolation' : 'second_order',\
        'wfc_extrapolation' : 'second_order'\
    }
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_nvt.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_npt(cell,xc,pot,tempw,tolp,press,dt=0.5,level='low',estep=1000,nstep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
            f = open(fpot,'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                if 'Suggested minimum cutoff for wavefunctions:' in line:
                    line = line.split()
                    if ecutwfc_ <= float(line[5]):
                        ecutwfc_ = float(line[5])
                if 'Suggested minimum cutoff for charge density:' in line:
                    line = line.split()
                    if ecutrho_ <= float(line[6]):
                        ecutrho_ = float(line[6])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'vc-md',\
        'dt'               : dt*1.0e-15/4.8378e-17,\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : True,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nstep'            : nstep,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'nosym'            : True,\
        'cell_dofree'      : 'all',\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'ecutrho'          : ecutrho_,\
        'ion_dynamics'     : 'beeman',\
        'ion_temperature'  : 'rescaling',\
        'tempw'            : tempw,\
        'tolp'             : tolp,\
        'cell_dynamics'    : 'w',\
        'press'            : press,\
        'pot_extrapolation' : 'second_order',\
        'wfc_extrapolation' : 'second_order'\
    }
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_npt.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_phonon(cell,nq=(4,4,4),nk=(16,16,16)):
    f = open('qe_ph.pwi','w')
    f.write('&INPUTPH'+'\n')
    f.write('outdir = "./outdir"'+'\n')
    f.write('tr2_ph =  1.0e-12'+'\n')
    f.write('fildyn  =  "ph.dyn"'+'\n')
    f.write('ldisp  =  .true.'+'\n')
    f.write('nq1  =  '+str(nq[0])+'\n')
    f.write('nq2  =  '+str(nq[1])+'\n')
    f.write('nq3  =  '+str(nq[2])+'\n')
    f.write('/'+'\n')
    f.close()
    f = open('qe_q2r.pwi','w')
    f.write('&INPUT'+'\n')
    f.write('fildyn  =  "ph.dyn"'+'\n')
    f.write('flfrc  =  "ph.fc"'+'\n')
    f.write('/'+'\n')
    f.close()
    f = open('qe_matdyn_phdos.pwi','w')
    f.write('&INPUT'+'\n')
    f.write('flfrc  =  "ph.fc"'+'\n')
    f.write('asr  =  "crystal"'+'\n')
    f.write('dos  =  .true.'+'\n')
    f.write('nk1  =  '+str(nk[0])+'\n')
    f.write('nk2  =  '+str(nk[1])+'\n')
    f.write('nk3  =  '+str(nk[2])+'\n')
    f.write('/'+'\n')
    f.close()
    f = open('qe_matdyn_phdisp.pwi','w')
    f.write('&INPUT'+'\n')
    f.write('flfrc  =  "ph.fc"'+'\n')
    f.write('asr  =  "crystal"'+'\n')
    f.write('q_in_band_form  =  .true.'+'\n')
    f.write('/'+'\n')
    lat = cell.cell.get_bravais_lattice()
    points = lat.get_special_points()
    paths = lat.special_path
    npath = len(paths)
    npoints = len(paths) - paths.count(',')
    f.write(str(npoints)+'\n')
    for n in range(npath):
        point = 20
        if paths[n] == ',':
            continue
        elif n+1 >= npath:
            point = 0
        elif n+1 < npath and paths[n+1] == ',':
            point = 0
        f.write(str(points[paths[n]][0])+' '+str(points[paths[n]][1])+' '+str(points[paths[n]][2])+' '+str(point)+' !'+paths[n]+'\n')
    f.close()

def set_qepot(cell,xc,pot0,level):
    from AtomicVirtuaLab.io import cell2atomlist
    if pot0 == 'us':
        pot = 'rrkjus'
    elif pot0 == 'paw':
        pot = 'kjpaw'
    if level == 'low':
        pseudo = \
            {'H' : 'H.'+str(xc)+'-'+str(pot)+'_psl.1.0.0.UPF',\
             'He': 'He.'+str(xc)+'-'+str(pot)+'_psl.1.0.0.UPF',\
             'Li': 'Li.'+str(xc)+'-sl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Be': 'Be.'+str(xc)+'-sl-'+str(pot)+'_psl.1.0.0.UPF',\
             'B' : 'B.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'C' : 'C.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'N' : 'N.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'O' : 'O.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'F' : 'F.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ne': 'Ne.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Na': 'Na.'+str(xc)+'-spnl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mg': 'Mg.'+str(xc)+'-spnl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Al': 'Al.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Si': 'Si.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'P' : 'P.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'S' : 'S.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cl': 'Cl.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ar': 'Ar.'+str(xc)+'-nl-'+str(pot)+'_psl.1.0.0.UPF',\
             'K' : 'K.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ca': 'Ca.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sc': 'Sc.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ti': 'Ti.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'V' : 'V.'+str(xc)+'-spnl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cr': 'Cr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mn': 'Mn.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Fe': 'Fe.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Co': 'Co.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ni': 'Ni.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cu': 'Cu.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Zn': 'Zn.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ga': 'Ga.'+str(xc)+'-dnl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ge': 'Ge.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'As': 'As.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Se': 'Se.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Br': 'Br.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Kr': 'Kr.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rb': 'Rb.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sr': 'Sr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Y' : 'Y.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Zr': 'Zr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Nb': 'Nb.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mo': 'Mo.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tc': 'Tc.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ru': 'Ru.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rh': 'Rh.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pd': 'Pd.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ag': 'Ag.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cd': 'Cd.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'In': 'In.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sn': 'Sn.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sb': 'Sb.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Te': 'Te.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'I' : 'I.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Xe': 'Xe.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cs': 'Cs.'+str(xc)+'-spnl-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ba': 'Ba.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'La': 'La.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ce': 'Ce.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pr': 'Pr.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Nd': 'Nd.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pm': 'Pm.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sm': 'Sm.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Eu': 'Eu.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Gd': 'Gd.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tb': 'Tb.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Dy': 'Dy.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ho': 'Ho.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Er': 'Er.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tm': 'Tm.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Yb': 'Yb.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Lu': 'Lu.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Hf': 'Hf.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ta': 'Ta.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'W' : 'W.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Re': 'Re.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Os': 'Os.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ir': 'Ir.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pt': 'Pt.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Au': 'Au.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Hg': 'Hg.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tl': 'Tl.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pb': 'Pb.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Bi': 'Bi.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Po': 'Po.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'At': 'At.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rn': 'Rn.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Fr': 'Fr.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ra': 'Ra.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ac': 'Ac.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Th': 'Th.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pa': 'Pa.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'U' : 'U.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Np': 'Np.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pu': 'Pu.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF'}
    elif level == 'high':
        pseudo = \
            {'H' : 'H.'+str(xc)+'-'+str(pot)+'_psl.1.0.0.UPF',\
             'He': 'He.'+str(xc)+'-'+str(pot)+'_psl.1.0.0.UPF',\
             'Li': 'Li.'+str(xc)+'-s-'+str(pot)+'_psl.1.0.0.UPF',\
             'Be': 'Be.'+str(xc)+'-s-'+str(pot)+'_psl.1.0.0.UPF',\
             'B' : 'B.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'C' : 'C.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'N' : 'N.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'O' : 'O.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'F' : 'F.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ne': 'Ne.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Na': 'Na.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mg': 'Mg.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Al': 'Al.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Si': 'Si.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'P' : 'P.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'S' : 'S.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cl': 'Cl.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ar': 'Ar.'+str(xc)+'-n-'+str(pot)+'_psl.1.0.0.UPF',\
             'K' : 'K.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ca': 'Ca.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sc': 'Sc.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ti': 'Ti.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'V' : 'V.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cr': 'Cr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mn': 'Mn.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Fe': 'Fe.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Co': 'Co.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ni': 'Ni.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cu': 'Cu.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Zn': 'Zn.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ga': 'Ga.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ge': 'Ge.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'As': 'As.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Se': 'Se.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Br': 'Br.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Kr': 'Kr.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rb': 'Rb.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sr': 'Sr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Y' : 'Y.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Zr': 'Zr.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Nb': 'Nb.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Mo': 'Mo.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tc': 'Tc.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ru': 'Ru.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rh': 'Rh.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pd': 'Pd.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ag': 'Ag.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cd': 'Cd.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'In': 'In.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sn': 'Sn.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sb': 'Sb.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Te': 'Te.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'I' : 'I.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Xe': 'Xe.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Cs': 'Cs.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ba': 'Ba.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'La': 'La.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ce': 'Ce.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pr': 'Pr.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Nd': 'Nd.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pm': 'Pm.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Sm': 'Sm.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Eu': 'Eu.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Gd': 'Gd.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tb': 'Tb.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Dy': 'Dy.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ho': 'Ho.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Er': 'Er.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tm': 'Tm.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Yb': 'Yb.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Lu': 'Lu.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Hf': 'Hf.'+str(xc)+'-spdfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ta': 'Ta.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'W' : 'W.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Re': 'Re.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Os': 'Os.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ir': 'Ir.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pt': 'Pt.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Au': 'Au.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Hg': 'Hg.'+str(xc)+'-spn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Tl': 'Tl.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pb': 'Pb.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Bi': 'Bi.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Po': 'Po.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'At': 'At.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Rn': 'Rn.'+str(xc)+'-dn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Fr': 'Fr.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ra': 'Ra.'+str(xc)+'-spdn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Ac': 'Ac.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Th': 'Th.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pa': 'Pa.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'U' : 'U.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Np': 'Np.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF',\
             'Pu': 'Pu.'+str(xc)+'-spfn-'+str(pot)+'_psl.1.0.0.UPF'}
    symbols = cell2atomlist(cell)
    pseudo0 = {}
    for symbol in symbols:
        pseudo0[symbol] = pseudo[symbol]
    return pseudo0

def get_valence(cell,pseudo,level):
    import AtomicVirtuaLab.globalv as g
    z_valence={}
    for symbol in pseudo:
        fpot = g.qepot+'/'+str(level)+'/'+str(pseudo[symbol])
        f = open(fpot,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            if 'z valence' in line.lower():
                z_valence[symbol] = float(line.split()[0])
            elif 'z_valence' in line.lower():
                if line.split()[0] == '<PP_HEADER':
                    line0 = list(filter(lambda x: 'z_valence' in x,
                                       line.split(' ')))[0]
                    z_valence[symbol] = float(line0.split('=')[-1].strip().strip('"'))
    valence=0
    magmom=[]
    for atom in cell:
        magmom.append(float(z_valence[atom.symbol]))
        valence = valence + int(z_valence[atom.symbol])
    return valence, magmom
