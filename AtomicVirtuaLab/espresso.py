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
    else:
        valence, magmom = get_valence(cell,pseudo,level)
        nbnd = int(valence/2.0)
        print('valence/2 = ', nbnd)
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

def plot_qe_band(result_dir,efermi,nspin=False):
    import re
    import sys
    from matplotlib import pyplot as plt
    
    if not nspin:
        f = open(result_dir+'/bands.out')
        l_bands_up = f.readlines()
        f.close()
    elif nspin:
        f = open(result_dir+'/bands_up.out')
        l_bands_up = f.readlines()
        f.close()
        f = open(result_dir+'/bands_down.out')
        l_bands_down = f.readlines()
        f.close()
    
    f = open(result_dir+'/qe_nscf_band.pwi','r')
    l_pwi = f.readlines()
    f.close()

    spoint=[]
    flg = False
    for line in l_pwi:
        if 'K_POINTS' in line:
            flg = True
            continue
        line = line.split()
        if flg and len(line) == 1:
            continue
        elif flg and len(line) == 0:
            flg = False
            continue
        elif flg:
            sp = line[4].replace('!', '')
            if sp == 'G':
                sp = u'\u0393'
            np = int(line[3])
            spoint.append([sp,np])

    tmp = l_bands_up[0].split(',')
    nbnd_up = int(re.sub(r'\D', '', tmp[0]))
    nks_up = int(re.sub(r'\D', '', tmp[1]))
    ibnd = 0
    iks = 0
    fks = True
    fband = False
    bands_up=[]
    for line in l_bands_up:
        if '&' in line:
            tmp1 = []
            continue
        line = line.split()
        if fks and not fband:
            tmp1.append([float(line[0]),float(line[1]),float(line[2])])
            fks = False
            fband = True
            tmp2 = []
        elif not fks and fband:
            for l in line:
                tmp2.append(float(l))
                ibnd = ibnd + 1
                if ibnd == nbnd_up:
                    bands_up.append([tmp1,tmp2])
                    ibnd = 0
                    fks = True
                    fband = False
                    tmp1 = []
    if nspin:
        tmp = l_bands_down[0].split(',')
        nbnd_down = int(re.sub(r'\D', '', tmp[0]))
        nks_down = int(re.sub(r'\D', '', tmp[1]))
        ibnd = 0
        iks = 0
        fks = True
        fband = False
        bands_down=[]
        for line in l_bands_down:
            if '&' in line:
                tmp1 = []
                continue
            line = line.split()
            if fks and not fband:
                tmp1.append([float(line[0]),float(line[1]),float(line[2])])
                fks = False
                fband = True
                tmp2 = []
            elif not fks and fband:
                for l in line:
                    tmp2.append(float(l))
                    ibnd = ibnd + 1
                    if ibnd == nbnd_down:
                        bands_down.append([tmp1,tmp2])
                        ibnd = 0
                        fks = True
                        fband = False
                        tmp1 = []

    x = []
    for n in range(nks_up):
        x.append(int(n))

    nsp = len(spoint)
    xtic = []
    xlab = []
    isp = 0
    j = 0
    distx=[]
    for i in x:
        if i == 0:
            xtic.append(int(i))
            xlab.append(spoint[i][0])
            j = i
            isp = isp + 1
            distx.append(i)
        elif i == j + spoint[isp-1][1] and spoint[isp][1] != 0:
            xtic.append(int(i))
            xlab.append(spoint[isp][0])
            j = i
            isp = isp + 1
        elif i == j + spoint[isp-1][1] and spoint[isp][1] == 0 and isp != nsp-1:
            xtic.append(int(i)+1)
            lab = spoint[isp][0]+'|'+spoint[isp+1][0]
            xlab.append(lab)
            j = i + 1
            isp = isp + 2
            distx.append(i+1)
        elif i == j + spoint[isp-1][1] and spoint[isp][1] == 0 and isp == nsp-1:
            xtic.append(int(i))
            xlab.append(spoint[isp][0])
            distx.append(i)

    for i in range(nbnd_up):
        y=[]
        for bands in bands_up:
            e = [r for r in bands[1]]
            y.append(e[i]-efermi)
        if not nspin:
            for j in range(len(distx)-1):
                plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]], color="#1d50a2", linewidth = 0.6)
        elif nspin:
            if i == 0:
                for j in range(len(distx)-1):
                    if j == 0:
                        plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]],label = 'up_spin', color="#1d50a2", linewidth = 0.6)
                    else:
                        plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]], color="#1d50a2", linewidth = 0.6)
            else:
                for j in range(len(distx)-1):
                    plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]], color="#1d50a2", linewidth = 0.6)
    if nspin:
        for i in range(nbnd_down):
            y=[]
            for bands in bands_down:
                e = [r for r in bands[1]]
                y.append(e[i]-efermi)
            if not nspin:
                for j in range(len(distx)-1):
                    plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]],linestyle = "dashed", color="#1d50a2", linewidth = 0.6)
            elif nspin:
                if i == 0:
                    for j in range(len(distx)-1):
                        if j == 0:
                            plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]],label = 'down_spin',linestyle = "dashed", color="#1d50a2", linewidth = 0.6)
                        else:
                            plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]],linestyle = "dashed", color="#1d50a2", linewidth = 0.6)
                else:
                    for j in range(len(distx)-1):
                        plt.plot(x[distx[j]:distx[j+1]],y[distx[j]:distx[j+1]],linestyle = "dashed", color="#1d50a2", linewidth = 0.6)
    plt.xlim(0,nks_up-1)
    plt.ylim(-10.0,10.0)
    plt.axhline(y=0.0, xmin=0.0, xmax=1.0,ls="dashed", color="black", linewidth = 0.6)
    plt.yticks(fontsize=8)
    plt.xticks(xtic,xlab,fontsize=8)
    plt.xlabel('Wave vector', fontsize=10)
    plt.ylabel('E-Efermi(eV)', fontsize=10)
    if nspin:
        plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    for xv in xtic:
        if xv == min(xtic) or xv == max(xtic):
            continue
        plt.axvline(x=xv, ymin=0.0, ymax=1.0, color="black", linewidth = 0.6)
    plt.savefig('qe_band.png',dpi=500)
    #plt.show()
    plt.clf()

def plot_qe_dos(result_dir,efermi,nspin=False):
    import os
    import re
    import sys
    from matplotlib import pyplot as plt
    import pandas as pd

    symbols = ['Si','O']

    files = os.listdir(result_dir+'/')
    fpdos = []
    for f in files:
        if 'pdos_atm#' in f and 'wfc#' in f:
            fpdos.append(f)
    d_pdos={}
    for symbol in symbols:
        l_pdos=[]
        for f in fpdos:
            tmp = f.split('_')
            for s in tmp:
                if 'atm#' in s and symbol in s:
                    l_pdos.append(f)
        d_pdos[symbol] = l_pdos
    pdos_up={}
    pdos_down={}
    ados_up={}
    ados_down={}
    for symbol in d_pdos:
        for f_pdos in d_pdos[symbol]:
            tmp = f_pdos.split('_')
            for s in tmp:
                if 'wfc#' in s:
                    wfc = s.replace('wfc#','')
                    wfc = re.split(r'[()]',wfc)
                    wfc.remove('')
                    if wfc[1] == 's':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1])
                        pdos_up[key]=[]
                        pdos_down[key]=[]
                    elif wfc[1] == 'p':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'y')
                        pdos_up[key] = []
                        pdos_down[key] = []
                    elif wfc[1] == 'd':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zx')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zy')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x2-y2')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'xy')
                        pdos_up[key] = []
                        pdos_down[key] = []
                    elif wfc[1] == 'f':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z3')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2x')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2y')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z(x2-y2)')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zxy')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x(x2-3y2)')
                        pdos_up[key] = []
                        pdos_down[key] = []
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'y(3x2-y2)')
                        pdos_up[key] = []
                        pdos_down[key] = []
    for symbol in d_pdos:
        for f_pdos in d_pdos[symbol]:
            f = open(result_dir+'/'+f_pdos,'r')
            lines = f.readlines()
            f.close()
            l_pdos = []
            for line in lines:
                if '#' in line:
                    continue
                line = line.split()
                l_pdos.append(line)
            E = [float(r[0])-efermi for r in l_pdos]
            tmp = f_pdos.split('_')
            for s in tmp:
                if 'wfc#' in s:
                    wfc = s.replace('wfc#','')
                    wfc = re.split(r'[()]',wfc)
                    wfc.remove('')
                    if wfc[1] == 's':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1])
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                    elif wfc[1] == 'p':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'y')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                    elif wfc[1] == 'd':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zx')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zy')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x2-y2')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[5]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[9]) for r in l_pdos]
                                p_down = [float(r[10]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[5]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[9]) for r in l_pdos]
                                p_down = [float(r[10]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'xy')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[6]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[11]) for r in l_pdos]
                                p_down = [float(r[12]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[6]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[11]) for r in l_pdos]
                                p_down = [float(r[12]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                    elif wfc[1] == 'f':
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z3')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[2]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[3]) for r in l_pdos]
                                p_down = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2x')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[3]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[5]) for r in l_pdos]
                                p_down = [float(r[6]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z2y')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[4]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[7]) for r in l_pdos]
                                p_down = [float(r[8]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'z(x2-y2)')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[5]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[9]) for r in l_pdos]
                                p_down = [float(r[10]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[5]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[9]) for r in l_pdos]
                                p_down = [float(r[10]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'zxy')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[6]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[11]) for r in l_pdos]
                                p_down = [float(r[12]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[6]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[11]) for r in l_pdos]
                                p_down = [float(r[12]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'x(x2-3y2)')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[7]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[13]) for r in l_pdos]
                                p_down = [float(r[14]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[7]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[13]) for r in l_pdos]
                                p_down = [float(r[14]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]
                        key = symbol+'_'+str(wfc[0])+str(wfc[1]+'y(3x2-y2)')
                        if len(pdos_up[key]) == 0:
                            if nspin == False:
                                p_up = [float(r[8]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                            elif nspin == True:
                                p_up = [float(r[15]) for r in l_pdos]
                                p_down = [float(r[16]) for r in l_pdos]
                                pdos_up[key] = [E,p_up].copy()
                                pdos_down[key] = [E,p_down].copy()
                        elif len(pdos_up[key]) != 0:
                            if nspin == False:
                                p_up = [float(r[8]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                            elif nspin == True:
                                p_up = [float(r[15]) for r in l_pdos]
                                p_down = [float(r[16]) for r in l_pdos]
                                pdos_up[key][1] = [x+y for x,y in zip(pdos_up[key][1],p_up)]
                                pdos_down[key][1] = [x+y for x,y in zip(pdos_down[key][1],p_down)]

    sum_pdos_up={}
    sum_pdos_down={}
    sum_ados_up={}
    sum_ados_down={}
    tot_dos_up={}
    tot_dos_down={}
    for symbol in symbols:
        sum_ados_up[symbol]=[]
        sum_ados_down[symbol]=[]
        for orb in ['s','pz','px','py','dz2','dzx','dzy','dx2-y2','dxy','fz3','fz2x','fz2y','fz(x2-y2)','fzxy','fx(x2-3y2)','fy(3x2-y2)']:
            key = symbol+'_'+orb
            sum_pdos_up[key] = []
            sum_pdos_down[key] = []
    tot_dos_up['tot'] = []
    tot_dos_down['tot'] = []

    for symbol in symbols:
        for orb in ['s','pz','px','py','dz2','dzx','dzy','dx2-y2','dxy','fz3','fz2x','fz2y','fz(x2-y2)','fzxy','fx(x2-3y2)','fy(3x2-y2)']:
            for pdos in pdos_up:
                pdos0 = pdos.split('_')
                if symbol == pdos0[0] and orb in pdos0[1]:
                    key = symbol+'_'+orb
                    if len(sum_pdos_up[key]) == 0:
                        sum_pdos_up[key] = pdos_up[pdos].copy()
                    elif len(sum_pdos_up[key]) != 0:
                        sum_pdos_up[key][1] = [x+y for x,y in zip(sum_pdos_up[key][1],pdos_up[pdos][1])]
                    if len(sum_ados_up[symbol]) == 0:
                        sum_ados_up[symbol] = pdos_up[pdos].copy()
                    elif len(sum_ados_up[symbol]) != 0:
                        sum_ados_up[symbol][1] = [x+y for x,y in zip(sum_ados_up[symbol][1],pdos_up[pdos][1])]
                    if len(tot_dos_up['tot']) == 0:
                        tot_dos_up['tot'] = pdos_up[pdos].copy()
                    elif len(tot_dos_up['tot']) != 0:
                        tot_dos_up['tot'][1] = [x+y for x,y in zip(tot_dos_up['tot'][1],pdos_up[pdos][1])]
            if nspin == True:
                for pdos in pdos_down:
                    pdos0 = pdos.split('_')
                    if symbol == pdos0[0] and orb in pdos0[1]:
                        key = symbol+'_'+orb
                        if len(sum_pdos_down[key]) == 0:
                            sum_pdos_down[key] = pdos_down[pdos].copy()
                        elif len(sum_pdos_down[key]) != 0:
                            sum_pdos_down[key][1] = [x+y for x,y in zip(sum_pdos_down[key][1],pdos_down[pdos][1])]
                        if len(sum_ados_down[symbol]) == 0:
                            sum_ados_down[symbol] = pdos_down[pdos].copy()
                        elif len(sum_ados_down[symbol]) != 0:
                            sum_ados_down[symbol][1] = [x+y for x,y in zip(sum_ados_down[symbol][1],pdos_down[pdos][1])]
                        if len(tot_dos_down['tot']) == 0:
                            tot_dos_down['tot'] = pdos_down[pdos].copy()
                        elif len(tot_dos_down['tot']) != 0:
                            tot_dos_down['tot'][1] = [x+y for x,y in zip(tot_dos_down['tot'][1],pdos_down[pdos][1])]
    
    # orbital DOS
    data = pd.DataFrame(pdos_up,index=('energy','pdos'))
    print(data.columns.values)
    for orb in data.columns.values:
        data2 = pd.DataFrame(data[orb]['pdos'],index=data[orb]['energy'],columns=[orb])
        data2.to_csv(orb+'.csv')
    if nspin == False:
        for pdos in pdos_up:
            plt.plot(pdos_up[pdos][0],pdos_up[pdos][1],label=pdos+'(UP)', linewidth = 0.6)
    if nspin == True:
        for pdos in pdos_up:
            plt.plot(pdos_up[pdos][0],pdos_up[pdos][1],label=pdos+'(UP)', linewidth = 0.6)
        for pdos in pdos_down:
            pdos_down[pdos][1] = [-1.0*x for x in pdos_down[pdos][1]]
            plt.plot(pdos_down[pdos][0],pdos_down[pdos][1],label=pdos+'(DOWN)', linewidth = 0.6)
    plt.xlabel('E-Efermi(eV)')
    plt.ylabel('DOS')
    #plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    plt.legend(loc='upper right',fontsize=8)
    plt.ylim(0.0,20.0)
    plt.xlim(-10.0,50.0)
    if nspin == False:
        plt.gca().set_ylim(bottom=0)
    plt.tick_params(labelsize=10)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black")
    plt.savefig('qe_orbital_dos.png',dpi=500)
    #plt.show()
    plt.clf()

    """
    #from matplotlib import pyplot as plt
    if nspin == False:
        plt.plot(tot_dos_up['tot'][0],tot_dos_up['tot'][1],label='total')
        for pdos in sum_pdos_up:
            if len(sum_pdos_up[pdos]) == 0:
                continue
            plt.plot(sum_pdos_up[pdos][0],sum_pdos_up[pdos][1],label=pdos)
    if nspin == True:
        plt.plot(tot_dos_up['tot'][0],tot_dos_up['tot'][1],label='total_up')
        tot_dos_down['tot'][1] = [-1.0*x for x in tot_dos_down['tot'][1]]
        plt.plot(tot_dos_down['tot'][0],tot_dos_down['tot'][1],label='total_down')
        for pdos in sum_pdos_up:
            if len(sum_pdos_up[pdos]) == 0:
                continue
            plt.plot(sum_pdos_up[pdos][0],sum_pdos_up[pdos][1],label=pdos+'_up')
        for pdos in sum_pdos_down:
            if len(sum_pdos_down[pdos]) == 0:
                continue
            sum_pdos_down[pdos][1] = [-1.0*x for x in sum_pdos_down[pdos][1]]
            plt.plot(sum_pdos_down[pdos][0],sum_pdos_down[pdos][1],label=pdos+'_down')
    plt.xlabel('E-Efermi(eV)')
    plt.ylabel('DOS')
    plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    plt.xlim(-10.0,10.0)
    if nspin == False:
        plt.gca().set_ylim(bottom=0)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black")
    plt.show()
    """

    #Total DOS
    plt.clf()
    if nspin == False:
        plt.plot(tot_dos_up['tot'][0],tot_dos_up['tot'][1],label='Total DOS(UP)', linewidth = 0.6)
        #for pdos in sum_ados_up:
        #    if len(sum_ados_up[pdos]) == 0:
        #        continue
        #    plt.plot(sum_ados_up[pdos][0],sum_ados_up[pdos][1],label=pdos)
    if nspin == True:
        plt.plot(tot_dos_up['tot'][0],tot_dos_up['tot'][1],label='Total DOS(UP)', linewidth = 0.6)
        #tot_dos_down['tot'][1] = [-1.0*x for x in tot_dos_down['tot'][1]]
        plt.plot(tot_dos_down['tot'][0],tot_dos_down['tot'][1],label='Total DOS(DOWN)', linewidth = 0.6)
        #for pdos in sum_ados_up:
        #    if len(sum_ados_up[pdos]) == 0:
        #        continue
        #    plt.plot(sum_ados_up[pdos][0],sum_ados_up[pdos][1],label=pdos+'_up')
        #for pdos in sum_ados_down:
        #    if len(sum_ados_down[pdos]) == 0:
        #        continue
        #    sum_ados_down[pdos][1] = [-1.0*x for x in sum_ados_down[pdos][1]]
        #    plt.plot(sum_ados_down[pdos][0],sum_ados_down[pdos][1],label=pdos+'_down')
    plt.xlabel('E-Efermi(eV)')
    plt.ylabel('DOS')
    #plt.legend(loc='upper left',bbox_to_anchor=(1,1))
    plt.legend(loc='upper right',fontsize=8)
    plt.ylim(0.0,20.0)
    plt.xlim(-5.0,5.0)
    if nspin == False:
        plt.gca().set_ylim(bottom=0)
    plt.tick_params(labelsize=10)
    plt.axvline(x=0.0, ymin=0.0, ymax=1.0,ls="dashed", color="black", linewidth = 0.6)
    plt.savefig('qe_total_dos.png',dpi=500)
    #plt.show()
    
    plt.clf()

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
