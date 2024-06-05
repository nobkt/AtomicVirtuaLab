def mk_qe_input_scf(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',tstress=True,nosym=False,options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'scf',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : tstress,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nosym'            : nosym,\
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
        nbnd = int(valence/2.0)
        #print('valence/2 = ', nbnd)
    cell.write('qe_scf.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe2yambo_input_scf(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',tstress=True,nosym=False,options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'scf',\
        'restart_mode'     : 'from_scratch',\
        'verbosity'        : 'high',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : tstress,\
        'tprnfor'          : False,\
        'outdir'           : './outdir',\
        'nosym'            : nosym,\
        'force_symmorphic' : True,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_\
    }
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
        nbnd = int(valence/2.0)
        print('valence/2 = ', nbnd)
    cell.write('qe_scf.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,koffset=(1,1,1),crystal_coordinates=False)

def mk_qe2yambo_input_scf_gaupbe(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',tstress=True,nosym=False,options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'scf',\
        'restart_mode'     : 'from_scratch',\
        'verbosity'        : 'high',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : tstress,\
        'tprnfor'          : False,\
        'outdir'           : './outdir',\
        'nosym'            : nosym,\
        'force_symmorphic' : True,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_,\
        'input_dft'        : 'gaupbe',\
        'x_gamma_extrapolation' : False,\
        'exxdiv_treatment' : 'none'\
    }
    if kpts == None:
        input_data['nqx1'] = 1
        input_data['nqx2'] = 1
        input_data['nqx3'] = 1
    else:
        input_data['nqx1'] = int(kpts[0])
        input_data['nqx2'] = int(kpts[1])
        input_data['nqx3'] = int(kpts[2])       
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
        nbnd = int(valence/2.0)
        print('valence/2 = ', nbnd)
    cell.write('qe_scf.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,koffset=(1,1,1),crystal_coordinates=False)


def mk_qe2yambo_input_nscf(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',tstress=True,nosym=False,options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'nscf',\
        'restart_mode'     : 'from_scratch',\
        'verbosity'        : 'high',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : tstress,\
        'tprnfor'          : True,\
        'outdir'           : './outdir',\
        'nosym'            : nosym,\
        'force_symmorphic' : True,\
        'occupations'      : 'smearing',\
        'smearing'         : 'gaussian',\
        'degauss'          : 0.01,\
        'conv_thr'         : 1.0e-6,\
        'diagonalization'  : 'david',\
        'mixing_beta'      : mixing_beta,\
        'electron_maxstep' : estep,\
        'ecutwfc'          : ecutwfc_\
    }
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
        nbnd = int(valence/2.0)
        print('valence/2 = ', nbnd)
    cell.write('qe_nscf.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)


def mk_qe_input_dos(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
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

def mk_qe_input_band(cell,xc,pot,level='low',estep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import sys
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    else:
        valence, magmom = get_valence(cell,pseudo,level,xc)
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

def mk_qe_input_relax(cell,xc,pot,level='low',estep=1000,nstep=1000,nosym=False,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',tstress=True,options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
    input_data={\
        'calculation'      : 'relax',\
        'restart_mode'     : 'from_scratch',\
        'wf_collect'       : True,\
        'pseudo_dir'       : './pseudo',\
        'tstress'          : tstress,\
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_relax.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_vcrelax(cell,xc,pot,level='low',estep=1000,nstep=1000,nosym=False,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_vc-relax.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_nvt(cell,xc,pot,tempw,tolp,dt=0.5,level='low',estep=1000,nstep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',options={},nspin=False):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
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
    if len(options) != 0:
        for option in options:
            input_data[option] = options[option]
    if nspin == True:
        input_data['nspin'] = 2
        valence, magmom = get_valence(cell,pseudo,level,xc)
        cell.set_initial_magnetic_moments(magmom)
        nbnd = int(valence/2.0*2.0)
        input_data['nbnd'] = nbnd
    cell.write('qe_nvt.pwi',input_data=input_data,pseudopotentials=pseudo,kpts=kpts,crystal_coordinates=False)

def mk_qe_input_npt(cell,xc,pot,tempw,tolp,press,dt=0.5,level='low',estep=1000,nstep=1000,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=None,ecut='manual',cell_dofree='all',options={},nspin=False,magElm=None):
    from ase.io import write
    import AtomicVirtuaLab.globalv as g
    import shutil
    import os
    import json
    level0 = level
    pseudo = set_qepot(cell,xc,pot,level0)
    if ecut == 'manual':
        ecutwfc_ = ecutwfc
        ecutrho_ = ecutrho
    elif ecut == 'auto' and level != 'SSSP_precision':
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
    elif ecut == 'auto' and level == 'SSSP_precision':
        ecutwfc_ = 0.0
        ecutrho_ = 0.0
        for symbol in pseudo:
            fjson = g.qepot+'/SSSP_1.3.0_PBE_precision.json'
            f = open(fjson,'r')
            json_load = json.load(f)
            f.close()
            if ecutwfc_ <= float(json_load[symbol]['cutoff_wfc']):
                ecutwfc_ = float(json_load[symbol]['cutoff_wfc'])
            if ecutrho_ <= float(json_load[symbol]['cutoff_rho']):
                ecutrho_ = float(json_load[symbol]['cutoff_rho'])
    os.makedirs('./pseudo',exist_ok=True)
    for symbol in pseudo:
        shutil.copy(g.qepot+'/'+str(level)+'_'+str(xc)+'/'+str(pseudo[symbol]),'./pseudo')
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
        'cell_dofree'      : str(cell_dofree),\
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
        valence, magmom0 = get_valence(cell,pseudo,level,xc)
        if magElm is not None:
            i = 0
            magmom=[]
            for atom in cell:
                if atom.symbol in magElm:
                    magmom.append(magmom0[i])
                else:
                    magmom.append(0)
                i = i + 1
            magmom0 = magmom    
        cell.set_initial_magnetic_moments(magmom0)
        #nbnd = int(valence/2.0*1.2)
        #input_data['nbnd'] = nbnd
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

def plot_qe_dos(result_dir,efermi,symbols,nspin=False):
    import os
    import re
    import sys
    from matplotlib import pyplot as plt
    import pandas as pd

    #symbols = ['Si','O']

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
    plt.xlim(-10.0,10.0)
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
    plt.xlim(-10.0,10.0)
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
    elif level == 'fhi':
        pseudo = \
            {'Ac' : 'Ac.pbe-mt_fhi.UPF',\
             'Db' : 'Db.pbe-mt_fhi.UPF',\
             'Mo' : 'Mo.pbe-mt_fhi.UPF',\
             'Sb' : 'Sb.pbe-mt_fhi.UPF',\
             'Ag' : 'Ag.pbe-mt_fhi.UPF',\
             'Er' : 'Er.pbe-mt_fhi.UPF',\
             'Mt' : 'Mt.pbe-mt_fhi.UPF',\
             'Sc' : 'Sc.pbe-mt_fhi.UPF',\
             'Al' : 'Al.pbe-mt_fhi.UPF',\
             'Es' : 'Es.pbe-mt_fhi.UPF',\
             'N'  : 'N.pbe-mt_fhi.UPF',\
             'Se' : 'Se.pbe-mt_fhi.UPF',\
             'Am' : 'Am.pbe-mt_fhi.UPF',\
             'F'  : 'F.pbe-mt_fhi.UPF',\
             'Na' : 'Na.pbe-mt_fhi.UPF',\
             'Sg' : 'Sg.pbe-mt_fhi.UPF',\
             'Ar' : 'Ar.pbe-mt_fhi.UPF',\
             'Fe' : 'Fe.pbe-mt_fhi.UPF',\
             'Nb' : 'Nb.pbe-mt_fhi.UPF',\
             'Si' : 'Si.pbe-mt_fhi.UPF',\
             'As' : 'As.pbe-mt_fhi.UPF',\
             'Fm' : 'Fm.pbe-mt_fhi.UPF',\
             'Nd' : 'Nd.pbe-mt_fhi.UPF',\
             'Sm' : 'Sm.pbe-mt_fhi.UPF',\
             'At' : 'At.pbe-mt_fhi.UPF',\
             'Fr' : 'Fr.pbe-mt_fhi.UPF',\
             'Ne' : 'Ne.pbe-mt_fhi.UPF',\
             'Sn' : 'Sn.pbe-mt_fhi.UPF',\
             'Au' : 'Au.pbe-mt_fhi.UPF',\
             'Ga' : 'Ga.pbe-mt_fhi.UPF',\
             'Ni' : 'Ni.pbe-mt_fhi.UPF',\
             'Sr' : 'Sr.pbe-mt_fhi.UPF',\
             'B'  : 'B.pbe-mt_fhi.UPF',\
             'Gd' : 'Gd.pbe-mt_fhi.UPF',\
             'Np' : 'Np.pbe-mt_fhi.UPF',\
             'Ta' : 'Ta.pbe-mt_fhi.UPF',\
             'Ba' : 'Ba.pbe-mt_fhi.UPF',\
             'Ge' : 'Ge.pbe-mt_fhi.UPF',\
             'O'  : 'O.pbe-mt_fhi.UPF',\
             'Tc' : 'Tc.pbe-mt_fhi.UPF',\
             'Be' : 'Be.pbe-mt_fhi.UPF',\
             'H'  : 'H.pbe-mt_fhi.UPF',\
             'Os' : 'Os.pbe-mt_fhi.UPF',\
             'Te' : 'Te.pbe-mt_fhi.UPF',\
             'Bh' : 'Bh.pbe-mt_fhi.UPF',\
             'He' : 'He.pbe-mt_fhi.UPF',\
             'P'  : 'P.pbe-mt_fhi.UPF',\
             'Ti' : 'Ti.pbe-mt_fhi.UPF',\
             'Bi' : 'Bi.pbe-mt_fhi.UPF',\
             'Hf' : 'Hf.pbe-mt_fhi.UPF',\
             'Pa' : 'Pa.pbe-mt_fhi.UPF',\
             'Tl' : 'Tl.pbe-mt_fhi.UPF',\
             'Bk' : 'Bk.pbe-mt_fhi.UPF',\
             'Hg' : 'Hg.pbe-mt_fhi.UPF',\
             'Pb' : 'Pb.pbe-mt_fhi.UPF',\
             'Tm' : 'Tm.pbe-mt_fhi.UPF',\
             'Br' : 'Br.pbe-mt_fhi.UPF',\
             'Hs' : 'Hs.pbe-mt_fhi.UPF',\
             'Pd' : 'Pd.pbe-mt_fhi.UPF',\
             'U'  : 'U.pbe-mt_fhi.UPF',\
             'C'  : 'C.pbe-mt_fhi.UPF',\
             'I'  : 'I.pbe-mt_fhi.UPF',\
             'Pm' : 'Pm.pbe-mt_fhi.UPF',\
             'V'  : 'V.pbe-mt_fhi.UPF',\
             'Ca' : 'Ca.pbe-mt_fhi.UPF',\
             'In' : 'In.pbe-mt_fhi.UPF',\
             'Po' : 'Po.pbe-mt_fhi.UPF',\
             'W'  : 'W.pbe-mt_fhi.UPF',\
             'Cd' : 'Cd.pbe-mt_fhi.UPF',\
             'Ir' : 'Ir.pbe-mt_fhi.UPF',\
             'Pt' : 'Pt.pbe-mt_fhi.UPF',\
             'Xe' : 'Xe.pbe-mt_fhi.UPF',\
             'Ce' : 'Ce.pbe-mt_fhi.UPF',\
             'K'  : 'K.pbe-mt_fhi.UPF',\
             'Ra' : 'Ra.pbe-mt_fhi.UPF',\
             'Y'  : 'Y.pbe-mt_fhi.UPF',\
             'Cf' : 'Cf.pbe-mt_fhi.UPF',\
             'Kr' : 'Kr.pbe-mt_fhi.UPF',\
             'Rb' : 'Rb.pbe-mt_fhi.UPF',\
             'Yb' : 'Yb.pbe-mt_fhi.UPF',\
             'Cl' : 'Cl.pbe-mt_fhi.UPF',\
             'Li' : 'Li.pbe-mt_fhi.UPF',\
             'Re' : 'Re.pbe-mt_fhi.UPF',\
             'Zn' : 'Zn.pbe-mt_fhi.UPF',\
             'Cm' : 'Cm.pbe-mt_fhi.UPF',\
             'Lr' : 'Lr.pbe-mt_fhi.UPF',\
             'Rf' : 'Rf.pbe-mt_fhi.UPF',\
             'Zr' : 'Zr.pbe-mt_fhi.UPF',\
             'Co' : 'Co.pbe-mt_fhi.UPF',\
             'Lu' : 'Lu.pbe-mt_fhi.UPF',\
             'Rh' : 'Rh.pbe-mt_fhi.UPF',\
             'Cr' : 'Cr.pbe-mt_fhi.UPF',\
             'Md' : 'Md.pbe-mt_fhi.UPF',\
             'Rn' : 'Rn.pbe-mt_fhi.UPF',\
             'Cs' : 'Cs.pbe-mt_fhi.UPF',\
             'Mg' : 'Mg.pbe-mt_fhi.UPF',\
             'Ru' : 'Ru.pbe-mt_fhi.UPF',\
             'Cu' : 'Cu.pbe-mt_fhi.UPF',\
             'Mn' : 'Mn.pbe-mt_fhi.UPF',\
             'S'  : 'S.pbe-mt_fhi.UPF'}
    elif level == 'SG15':
        pseudo = \
            {'Ag' : 'Ag_ONCV_PBE-1.2.upf',\
             'Al' : 'Al_ONCV_PBE-1.2.upf',\
             'Ar' : 'Ar_ONCV_PBE-1.2.upf',\
             'As' : 'As_ONCV_PBE-1.2.upf',\
             'Au' : 'Au_ONCV_PBE-1.2.upf',\
             'B'  : 'B_ONCV_PBE-1.2.upf',\
             'Ba' : 'Ba_ONCV_PBE-1.2.upf',\
             'Be' : 'Be_ONCV_PBE-1.2.upf',\
             'Bi' : 'Bi_ONCV_PBE-1.2.upf',\
             'Br' : 'Br_ONCV_PBE-1.2.upf',\
             'C'  : 'C_ONCV_PBE-1.2.upf',\
             'Ca' : 'Ca_ONCV_PBE-1.2.upf',\
             'Cd' : 'Cd_ONCV_PBE-1.2.upf',\
             'Cl' : 'Cl_ONCV_PBE-1.2.upf',\
             'Co' : 'Co_ONCV_PBE-1.2.upf',\
             'Cr' : 'Cr_ONCV_PBE-1.2.upf',\
             'Cs' : 'Cs_ONCV_PBE-1.2.upf',\
             'Cu' : 'Cu_ONCV_PBE-1.2.upf',\
             'F'  : 'F_ONCV_PBE-1.2.upf',\
             'Fe' : 'Fe_ONCV_PBE-1.2.upf',\
             'Ga' : 'Ga_ONCV_PBE-1.2.upf',\
             'Ge' : 'Ge_ONCV_PBE-1.2.upf',\
             'H'  : 'H_ONCV_PBE-1.2.upf',\
             'He' : 'He_ONCV_PBE-1.2.upf',\
             'Hf' : 'Hf_ONCV_PBE-1.2.upf',\
             'Hg' : 'Hg_ONCV_PBE-1.2.upf',\
             'I'  : 'I_ONCV_PBE-1.2.upf',\
             'In' : 'In_ONCV_PBE-1.2.upf',\
             'Ir' : 'Ir_ONCV_PBE-1.2.upf',\
             'K'  : 'K_ONCV_PBE-1.2.upf',\
             'Kr' : 'Kr_ONCV_PBE-1.2.upf',\
             'La' : 'La_ONCV_PBE-1.2.upf',\
             'Li' : 'Li_ONCV_PBE-1.2.upf',\
             'Mg' : 'Mg_ONCV_PBE-1.2.upf',\
             'Mn' : 'Mn_ONCV_PBE-1.2.upf',\
             'Mo' : 'Mo_ONCV_PBE-1.2.upf',\
             'N'  : 'N_ONCV_PBE-1.2.upf',\
             'Na' : 'Na_ONCV_PBE-1.2.upf',\
             'Nb' : 'Nb_ONCV_PBE-1.2.upf',\
             'Ne' : 'Ne_ONCV_PBE-1.2.upf',\
             'Ni' : 'Ni_ONCV_PBE-1.2.upf',\
             'O'  : 'O_ONCV_PBE-1.2.upf',\
             'Os' : 'Os_ONCV_PBE-1.2.upf',\
             'P'  : 'P_ONCV_PBE-1.2.upf',\
             'Pb' : 'Pb_ONCV_PBE-1.2.upf',\
             'Pd' : 'Pd_ONCV_PBE-1.2.upf',\
             'Pt' : 'Pt_ONCV_PBE-1.2.upf',\
             'Rb' : 'Rb_ONCV_PBE-1.2.upf',\
             'Re' : 'Re_ONCV_PBE-1.2.upf',\
             'Rh' : 'Rh_ONCV_PBE-1.2.upf',\
             'Ru' : 'Ru_ONCV_PBE-1.2.upf',\
             'S'  : 'S_ONCV_PBE-1.2.upf',\
             'Sb' : 'Sb_ONCV_PBE-1.2.upf',\
             'Sc' : 'Sc_ONCV_PBE-1.2.upf',\
             'Se' : 'Se_ONCV_PBE-1.2.upf',\
             'Si' : 'Si_ONCV_PBE-1.2.upf',\
             'Sn' : 'Sn_ONCV_PBE-1.2.upf',\
             'Sr' : 'Sr_ONCV_PBE-1.2.upf',\
             'Ta' : 'Ta_ONCV_PBE-1.2.upf',\
             'Tc' : 'Tc_ONCV_PBE-1.2.upf',\
             'Te' : 'Te_ONCV_PBE-1.2.upf',\
             'Ti' : 'Ti_ONCV_PBE-1.2.upf',\
             'Tl' : 'Tl_ONCV_PBE-1.2.upf',\
             'V'  : 'V_ONCV_PBE-1.2.upf',\
             'W'  : 'W_ONCV_PBE-1.2.upf',\
             'Xe' : 'Xe_ONCV_PBE-1.2.upf',\
             'Y'  : 'Y_ONCV_PBE-1.2.upf',\
             'Zn' : 'Zn_ONCV_PBE-1.2.upf',\
             'Zr' : 'Zr_ONCV_PBE-1.2.upf'}
    elif level == 'SSSP_precision':
        if xc == 'pbe':
            pseudo = \
                {'Ac' : 'Ac.us.z_11.ld1.psl.v1.0.0-high.upf',\
                 'Ag' : 'Ag_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Al' : 'Al.pbe-n-kjpaw_psl.1.0.0.UPF',\
                 'Am' : 'Am.paw.z_17.ld1.uni-marburg.v0.upf',\
                 'Ar' : 'Ar.paw.z_8.ld1.psl.v1.0.0-high.upf',\
                 'As' : 'As.nc.z_15.oncvpsp3.dojo.v4-std.upf',\
                 'At' : 'At.us.z_17.ld1.psl.v1.0.0-high.upf',\
                 'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf',\
                 'B'  : 'B_pbe_v1.01.uspp.F.UPF',\
                 'Ba' : 'Ba.nc.z_10.oncvpsp4.dojo.v4-sp.upf',\
                 'Be' : 'Be_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Bi' : 'Bi_pbe_v1.uspp.F.UPF',\
                 'Bk' : 'Bk.paw.z_19.ld1.uni-marburg.v0.upf',\
                 'C'  : 'C.pbe-n-kjpaw_psl.1.0.0.UPF',\
                 'Ca' : 'Ca_pbe_v1.uspp.F.UPF',\
                 'Cd' : 'Cd.paw.z_20.ld1.psl.v1.0.0-high.upf',\
                 'Ce' : 'Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf',\
                 'Cf' : 'Cf.paw.z_20.ld1.uni-marburg.v0.upf',\
                 'Cl' : 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF',\
                 'Cm' : 'Cm.paw.z_18.ld1.uni-marburg.v0.upf',\
                 'Co' : 'Co_pbe_v1.2.uspp.F.UPF',\
                 'Cs' : 'Cs.nc.z_9.oncvpsp3.dojo.v4-str.upf',\
                 'Cu' : 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf',\
                 'Dy' : 'Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf',\
                 'Er' : 'Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf',\
                 'Es' : 'Es.paw.z_21.ld1.uni-marburg.v0.upf',\
                 'Eu' : 'Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf',\
                 'F'  : 'F.oncvpsp.upf',\
                 'Fe' : 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF',\
                 'Fm' : 'Fm.paw.z_22.ld1.uni-marburg.v0.upf',\
                 'Fr' : 'Fr.paw.z_19.ld1.psl.v1.0.0-high.upf',\
                 'Ga' : 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF',\
                 'Gd' : 'Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf',\
                 'H'  : 'H_ONCV_PBE-1.0.oncvpsp.upf',\
                 'He' : 'He_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Hf' : 'Hf-sp.oncvpsp.upf',\
                 'Hg' : 'Hg.us.z_12.uspp.gbrv.v1.upf',\
                 'Ho' : 'Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf',\
                 'I'  : 'I.nc.z_17.oncvpsp4.sg15.v0.upf',\
                 'In' : 'In.pbe-dn-rrkjus_psl.0.2.2.UPF',\
                 'Ir' : 'Ir.us.z_31.ld1.psl.v1.0.0-high.upf',\
                 'K'  : 'K.pbe-spn-kjpaw_psl.1.0.0.UPF',\
                 'Kr' : 'Kr.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'La' : 'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf',\
                 'Lr' : 'Lr.paw.z_25.ld1.uni-marburg.v0.upf',\
                 'Lu' : 'Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf',\
                 'Md' : 'Md.paw.z_23.ld1.uni-marburg.v0.upf',\
                 'Mo' : 'Mo_ONCV_PBE-1.0.oncvpsp.upf',\
                 'N'  : 'N.oncvpsp.upf',\
                 'Na' : 'Na.paw.z_9.ld1.psl.v1.0.0-low.upf',\
                 'Nb' : 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',\
                 'Nd' : 'Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf',\
                 'Ne' : 'Ne.paw.z_8.ld1.psl.v1.0.0-high.upf',\
                 'No' : 'No.paw.z_24.ld1.uni-marburg.v0.upf',\
                 'Np' : 'Np.paw.z_15.ld1.uni-marburg.v0.upf',\
                 'O'  : 'O.pbe-n-kjpaw_psl.0.1.UPF',\
                 'Os' : 'Os_pbe_v1.2.uspp.F.UPF',\
                 'P'  : 'P.pbe-n-rrkjus_psl.1.0.0.UPF',\
                 'Pa' : 'Pa.paw.z_13.ld1.uni-marburg.v0.upf',\
                 'Pb' : 'Pb.pbe-dn-kjpaw_psl.0.2.2.UPF',\
                 'Pd' : 'Pd_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Pm' : 'Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf',\
                 'Po' : 'Po.pbe-dn-rrkjus_psl.1.0.0.UPF',\
                 'Pr' : 'Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf',\
                 'Pt' : 'Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF',\
                 'Pu' : 'Pu.paw.z_16.ld1.uni-marburg.v0.upf',\
                 'Ra' : 'Ra.paw.z_20.ld1.psl.v1.0.0-high.upf',\
                 'Rb' : 'Rb_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Re' : 'Re_pbe_v1.2.uspp.F.UPF',\
                 'Rh' : 'Rh_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Rn' : 'Rn.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'Ru' : 'Ru_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Sc' : 'Sc.pbe-spn-kjpaw_psl.0.2.3.UPF',\
                 'Se' : 'Se_pbe_v1.uspp.F.UPF',\
                 'Si' : 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',\
                 'Sm' : 'Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf',\
                 'Sn' : 'Sn_pbe_v1.uspp.F.UPF',\
                 'Sr' : 'Sr_pbe_v1.uspp.F.UPF',\
                 'Ta' : 'Ta_pbe_v1.uspp.F.UPF',\
                 'Tb' : 'Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf',\
                 'Tc' : 'Tc_ONCV_PBE-1.0.oncvpsp.upf',\
                 'Te' : 'Te.us.z_6.ld1.psl.v1.0.0-low.upf',\
                 'Th' : 'Th.paw.z_12.ld1.uni-marburg.v0.upf',\
                 'Tl' : 'Tl_pbe_v1.2.uspp.F.UPF',\
                 'Tm' : 'Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf',\
                 'U'  : 'U.paw.z_14.ld1.uni-marburg.v0.upf',\
                 'W'  : 'W_pbe_v1.2.uspp.F.UPF',\
                 'Xe' : 'Xe.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'Y'  : 'Y_pbe_v1.uspp.F.UPF',\
                 'Yb' : 'Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf',\
                 'Zn' : 'Zn_pbe_v1.uspp.F.UPF',\
                 'Zr' : 'Zr_pbe_v1.uspp.F.UPF',\
                 'Br' : 'br_pbe_v1.4.uspp.F.UPF',\
                 'Cr' : 'cr_pbe_v1.5.uspp.F.UPF',\
                 'Ge' : 'ge_pbe_v1.4.uspp.F.UPF',\
                 'Li' : 'li_pbe_v1.4.uspp.F.UPF',\
                 'Mg' : 'mg_pbe_v1.4.uspp.F.UPF',\
                 'Mn' : 'mn_pbe_v1.5.uspp.F.UPF',\
                 'Ni' : 'ni_pbe_v1.4.uspp.F.UPF',\
                 'S'  : 's_pbe_v1.4.uspp.F.UPF',\
                 'Sb' : 'sb_pbe_v1.4.uspp.F.UPF',\
                 'Ti' : 'ti_pbe_v1.4.uspp.F.UPF',\
                 'V'  : 'v_pbe_v1.4.uspp.F.UPF'}
        elif xc == 'pbesol':
            pseudo = \
                {'Ac' : 'Ac.us.z_11.ld1.psl.v1.0.0-high.upf',\
                 'Ag' : 'Ag_ONCV_PBEsol-1.0.upf',\
                 'Al' : 'Al.pbesol-n-kjpaw_psl.1.0.0.UPF',\
                 'Am' : 'Am.paw.z_17.ld1.uni-marburg.v0.upf',\
                 'Ar' : 'Ar.paw.z_8.ld1.psl.v1.0.0-high.upf',\
                 'As' : 'As.nc.z_15.oncvpsp3.dojo.v4-std.upf',\
                 'At' : 'At.us.z_17.ld1.psl.v1.0.0-high.upf',\
                 'Au' : 'Au_ONCV_PBEsol-1.0.upf',\
                 'Ba' : 'Ba.nc.z_10.oncvpsp4.dojo.v4-sp.upf',\
                 'Be' : 'Be_ONCV_PBEsol-1.0.upf',\
                 'Bk' : 'Bk.paw.z_19.ld1.uni-marburg.v0.upf',\
                 'C'  : 'C.pbesol-n-kjpaw_psl.1.0.0.UPF',\
                 'Cd' : 'Cd.paw.z_20.ld1.psl.v1.0.0-high.upf',\
                 'Ce' : 'Ce.paw.z_12.atompaw.wentzcovitch.v1.2.upf',\
                 'Cf' : 'Cf.paw.z_20.ld1.uni-marburg.v0.upf',\
                 'Cl' : 'Cl.pbesol-n-rrkjus_psl.1.0.0.UPF',\
                 'Cm' : 'Cm.paw.z_18.ld1.uni-marburg.v0.upf',\
                 'Cs' : 'Cs.nc.z_9.oncvpsp3.dojo.v4-str.upf',\
                 'Cu' : 'Cu.paw.z_11.ld1.psl.v1.0.0-low.upf',\
                 'Dy' : 'Dy.paw.z_20.atompaw.wentzcovitch.v1.2.upf',\
                 'Er' : 'Er.paw.z_22.atompaw.wentzcovitch.v1.2.upf',\
                 'Es' : 'Es.paw.z_21.ld1.uni-marburg.v0.upf',\
                 'Eu' : 'Eu.paw.z_17.atompaw.wentzcovitch.v1.2.upf',\
                 'F'  : 'F.oncvpsp.upf',\
                 'Fe' : 'Fe.pbesol-spn-kjpaw_psl.0.2.1.UPF',\
                 'Fm' : 'Fm.paw.z_22.ld1.uni-marburg.v0.upf',\
                 'Fr' : 'Fr.paw.z_19.ld1.psl.v1.0.0-high.upf',\
                 'Ga' : 'Ga.pbesol-dn-kjpaw_psl.1.0.0.UPF',\
                 'Gd' : 'Gd.paw.z_18.atompaw.wentzcovitch.v1.2.upf',\
                 'H'  : 'H_ONCV_PBEsol-1.0.upf',\
                 'He' : 'He_ONCV_PBEsol-1.0.upf',\
                 'Hf' : 'Hf-sp.oncvpsp.upf',\
                 'Hg' : 'Hg.us.z_12.uspp.gbrv.v1.upf',\
                 'Ho' : 'Ho.paw.z_21.atompaw.wentzcovitch.v1.2.upf',\
                 'I'  : 'I.nc.z_17.oncvpsp4.sg15.v0.upf',\
                 'In' : 'In.pbesol-dn-rrkjus_psl.0.2.2.UPF',\
                 'Ir' : 'Ir.us.z_31.ld1.psl.v1.0.0-high.upf',\
                 'K'  : 'K.pbesol-spn-kjpaw_psl.1.0.0.UPF',\
                 'Kr' : 'Kr.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'La' : 'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf',\
                 'Lr' : 'Lr.paw.z_25.ld1.uni-marburg.v0.upf',\
                 'Lu' : 'Lu.paw.z_25.atompaw.wentzcovitch.v1.2.upf',\
                 'Md' : 'Md.paw.z_23.ld1.uni-marburg.v0.upf',\
                 'Mo' : 'Mo_ONCV_PBEsol-1.0.upf',\
                 'N'  : 'N.oncvpsp.upf',\
                 'Na' : 'Na.paw.z_9.ld1.psl.v1.0.0-low.upf',\
                 'Nb' : 'Nb.pbesol-spn-kjpaw_psl.0.3.0.UPF',\
                 'Nd' : 'Nd.paw.z_14.atompaw.wentzcovitch.v1.2.upf',\
                 'Ne' : 'Ne.paw.z_8.ld1.psl.v1.0.0-high.upf',\
                 'No' : 'No.paw.z_24.ld1.uni-marburg.v0.upf',\
                 'Np' : 'Np.paw.z_15.ld1.uni-marburg.v0.upf',\
                 'O'  : 'O.pbesol-n-kjpaw_psl.0.1.UPF',\
                 'P'  : 'P.pbesol-n-rrkjus_psl.1.0.0.UPF',\
                 'Pa' : 'Pa.paw.z_13.ld1.uni-marburg.v0.upf',\
                 'Pb' : 'Pb.pbesol-dn-kjpaw_psl.0.2.2.UPF',\
                 'Pd' : 'Pd_ONCV_PBEsol-1.0.upf',\
                 'Pm' : 'Pm.paw.z_15.atompaw.wentzcovitch.v1.2.upf',\
                 'Po' : 'Po.pbesol-dn-rrkjus_psl.1.0.0.UPF',\
                 'Pr' : 'Pr.paw.z_13.atompaw.wentzcovitch.v1.2.upf',\
                 'Pt' : 'Pt.pbesol-spfn-rrkjus_psl.1.0.0.UPF',\
                 'Pu' : 'Pu.paw.z_16.ld1.uni-marburg.v0.upf',\
                 'Ra' : 'Ra.paw.z_20.ld1.psl.v1.0.0-high.upf',\
                 'Rb' : 'Rb_ONCV_PBEsol-1.0.upf',\
                 'Rh' : 'Rh_ONCV_PBEsol-1.0.upf',\
                 'Rn' : 'Rn.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'Ru' : 'Ru_ONCV_PBEsol-1.0.upf',\
                 'Sc' : 'Sc.pbesol-spn-kjpaw_psl.0.2.3.UPF',\
                 'Si' : 'Si.pbesol-n-rrkjus_psl.1.0.0.UPF',\
                 'Sm' : 'Sm.paw.z_16.atompaw.wentzcovitch.v1.2.upf',\
                 'Tb' : 'Tb.paw.z_19.atompaw.wentzcovitch.v1.2.upf',\
                 'Tc' : 'Tc_ONCV_PBEsol-1.0.upf',\
                 'Te' : 'Te.us.z_6.ld1.psl.v1.0.0-low.upf',\
                 'Th' : 'Th.paw.z_12.ld1.uni-marburg.v0.upf',\
                 'Tm' : 'Tm.paw.z_23.atompaw.wentzcovitch.v1.2.upf',\
                 'U'  : 'U.paw.z_14.ld1.uni-marburg.v0.upf',\
                 'Xe' : 'Xe.paw.z_18.ld1.psl.v1.0.0-high.upf',\
                 'Yb' : 'Yb.paw.z_24.atompaw.wentzcovitch.v1.2.upf',\
                 'B'  : 'b_pbesol_v1.4.uspp.F.UPF',\
                 'Bi' : 'bi_pbesol_v1.uspp.F.UPF',\
                 'Br' : 'br_pbesol_v1.4.uspp.F.UPF',\
                 'Ca' : 'ca_pbesol_v1.uspp.F.UPF',\
                 'Co' : 'co_pbesol_v1.2.uspp.F.UPF',\
                 'Cr' : 'cr_pbesol_v1.5.uspp.F.UPF',\
                 'Ge' : 'ge_pbesol_v1.4.uspp.F.UPF',\
                 'Li' : 'li_pbesol_v1.4.uspp.F.UPF',\
                 'Mg' : 'mg_pbesol_v1.4.uspp.F.UPF',\
                 'Mn' : 'mn_pbesol_v1.5.uspp.F.UPF',\
                 'Ni' : 'ni_pbesol_v1.4.uspp.F.UPF',\
                 'Os' : 'os_pbesol_v1.2.uspp.F.UPF',\
                 'Re' : 're_pbesol_v1.2.uspp.F.UPF',\
                 'S'  : 's_pbesol_v1.4.uspp.F.UPF',\
                 'Sb' : 'sb_pbesol_v1.4.uspp.F.UPF',\
                 'Se' : 'se_pbesol_v1.uspp.F.UPF',\
                 'Sn' : 'sn_pbesol_v1.4.uspp.F.UPF',\
                 'Sr' : 'sr_pbesol_v1.uspp.F.UPF',\
                 'Ta' : 'ta_pbesol_v1.uspp.F.UPF',\
                 'Ti' : 'ti_pbesol_v1.4.uspp.F.UPF',\
                 'Tl' : 'tl_pbesol_v1.2.uspp.F.UPF',\
                 'V'  : 'v_pbesol_v1.4.uspp.F.UPF',\
                 'W'  : 'w_pbesol_v1.2.uspp.F.UPF',\
                 'Y'  : 'y_pbesol_v1.4.uspp.F.UPF',\
                 'Zn' : 'zn_pbesol_v1.uspp.F.UPF',\
                 'Zr' : 'zr_pbesol_v1.uspp.F.UPF'}
    symbols = cell2atomlist(cell)
    pseudo0 = {}
    for symbol in symbols:
        pseudo0[symbol] = pseudo[symbol]
    return pseudo0

def get_valence(cell,pseudo,level,xc):
    import AtomicVirtuaLab.globalv as g
    if level == 'SSSP_precision':
        level = level+'_'+xc
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
                elif line.split('"')[0] == 'z_valence=':
                    line0 = line.split('"')
                    z_valence[symbol] = float(line0[1])
                else:
                    line0 = line.split('"')
                    z_valence[symbol] = float(line0[1])                    
    valence=0
    magmom=[]
    for atom in cell:
        magmom.append(float(z_valence[atom.symbol]))
        valence = valence + int(z_valence[atom.symbol])
    return valence, magmom
