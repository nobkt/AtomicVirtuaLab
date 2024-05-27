def mk_siesta_input_npt(cell,xc,basis_set,mesh_cutoff,kpts,temp,press,nstep,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={
                             'MaxSCFIterations':int(MaxSCFIterations),
                             'WriteCoorStep':True,
                             'WriteForces':True,
                             'WriteMDHistory':True,
                             'WriteCoorXmol':True,
                             'WriteCoorCerius':True,
                             'WriteMDXmol':True,
                             'MD.UseSaveXV':True,
                             'MD.TypeOfRun':'NoseParrinelloRahman',
                             'MD.ParrinelloRahmanMass':'1000 Ry*fs**2',
                             'MD.Target-pressure':str(press)+' GPa',
                             'MD.Target-temperature':str(temp)+' K',
                             'MD.Initial-temperature':str(temp)+' K',
                             'MD.Initial.Time.Step':1,
                             'MD.Final.Time.Step':int(nstep),
                             'MD.Length.Time.Step':str(dt)+' fs',
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':True,
                             #'WriteMullikenPop':1,
                             #'WriteHirshfeldPop':True,
                             #'WriteVoronoiPop':True,
                             #'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_nvt(cell,xc,basis_set,mesh_cutoff,kpts,temp,nstep,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,dt=0.5,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={
                             'MaxSCFIterations':int(MaxSCFIterations),
                             'WriteCoorStep':True,
                             'WriteForces':True,
                             'WriteMDHistory':True,
                             'WriteCoorXmol':True,
                             'WriteCoorCerius':True,
                             'WriteMDXmol':True,
                             'MD.UseSaveXV':True,
                             'MD.TypeOfRun':'Nose',
                             'MD.Target-temperature':str(temp)+' K',
                             'MD.Initial-temperature':str(temp)+' K',
                             'MD.Initial.Time.Step':1,
                             'MD.Final.Time.Step':int(nstep),
                             'MD.Length.Time.Step':str(dt)+' fs',
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':True,
                             #'WriteMullikenPop':1,
                             #'WriteHirshfeldPop':True,
                             #'WriteVoronoiPop':True,
                             #'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_optimize(cell,xc,basis_set,mesh_cutoff,kpts,nstep,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry

    fdf_arguments0={
                    'MaxSCFIterations':int(MaxSCFIterations),
                    'WriteCoorStep':True,
                    'WriteForces':True,
                    'WriteMDHistory':True,
                    'WriteCoorXmol':True,
                    'WriteCoorCerius':True,
                    'WriteMDXmol':True,
                    'MD.UseSaveXV':True,
                    'MD.TypeOfRun':'CG',
                    'MD.NumCGsteps':int(nstep),
                    'SolutionMethod':str(SolutionMethod),
                    'SCFMustConverge':True,
                    'WriteMullikenPop':1,
                    'WriteHirshfeldPop':True,
                    'WriteVoronoiPop':True,
                    'PartialChargesAtEveryGeometry':True,
                    'SCF.Mixers':['broyden'],
                    'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                    }
    
    if len(options) != 0:
        for option in options:
            fdf_arguments0[option] = options[option]

    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments=fdf_arguments0
              ).write_input(cell,'siesta')

def mk_siesta_input_cellopt(cell,xc,basis_set,mesh_cutoff,kpts,nstep,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={
                             'MaxSCFIterations':int(MaxSCFIterations),
                             'WriteCoorStep':True,
                             'WriteForces':True,
                             'WriteMDHistory':True,
                             'WriteCoorXmol':True,
                             'WriteCoorCerius':True,
                             'WriteMDXmol':True,
                             'MD.UseSaveXV':True,
                             'MD.TypeOfRun':'CG',
                             'MD.VariableCell':True,
                             'MD.NumCGsteps':int(nstep),
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':True,
                             #'WriteMullikenPop':1,
                             #'WriteHirshfeldPop':True,
                             #'WriteVoronoiPop':True,
                             #'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_scf(cell,xc,basis_set,mesh_cutoff,kpts,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,options={},spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    fdf_arguments0={
                   'MaxSCFIterations':int(MaxSCFIterations),
                    'WriteForces':True,
                    'WriteMDHistory':True,
                    'WriteCoorXmol':True,
                    'WriteCoorCerius':True,
                    'WriteMDXmol':True,
                    'MD.UseSaveXV':True,
                    'SolutionMethod':str(SolutionMethod),
                    'SCFMustConverge':True,
                    'WriteMullikenPop':1,
                    'WriteHirshfeldPop':True,
                    'WriteVoronoiPop':True,
                    'PartialChargesAtEveryGeometry':True,
                    'SCF.Mixers':['broyden'],
                    'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                    }
    
    if len(options) != 0:
        for option in options:
            fdf_arguments0[option] = options[option]

    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments=fdf_arguments0
              ).write_input(cell,'siesta')
    

def mk_siesta_input_scf_withEfield(cell,xc,basis_set,mesh_cutoff,kpts,pseudo_path,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={
                             'MaxSCFIterations':int(MaxSCFIterations),
                             'WriteCoorStep':True,
                             'WriteForces':True,
                             'WriteMDHistory':True,
                             'WriteCoorXmol':True,
                             'WriteCoorCerius':True,
                             'WriteMDXmol':True,
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':True,
                             #'WriteMullikenPop':1,
                             #'WriteHirshfeldPop':True,
                             #'WriteVoronoiPop':True,
                             #'PartialChargesAtEveryGeometry':True,
                             #'ExternalElectricField':[str(ex)+' '+str(ey)+' '+str(ez)+' V/Ang'],
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_scf_withEfield_wannier(cell,xc,basis_set,mesh_cutoff,kpts,pseudo_path,bandscale=1.0,ex=0.0,ey=0.0,ez=0.0,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    import os
    import sys
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={
                             'MaxSCFIterations':int(MaxSCFIterations),
                             'WriteCoorStep':True,
                             'WriteForces':True,
                             'WriteMDHistory':True,
                             'WriteCoorXmol':True,
                             'WriteCoorCerius':True,
                             'WriteMDXmol':True,
                             'SolutionMethod':str(SolutionMethod),
                             'SCF.DM.Tolerance':1E-08,
                             'SCFMustConverge':True,
                             #'WriteMullikenPop':1,
                             #'WriteHirshfeldPop':True,
                             #'WriteVoronoiPop':True,
                             #'PartialChargesAtEveryGeometry':True,
                             'ExternalElectricField':[str(ex)+' '+str(ey)+' '+str(ez)+' V/Ang'],
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)],
                             'siesta2Wannier90.WriteMmn':True,
                             'siesta2Wannier90.WriteAmn':True,
                             'siesta2Wannier90.WriteEig':True,
                             'siesta2Wannier90.WriteUnk':False
                             }
              ).write_input(cell,'siesta')
    
    dnelc = {}
    for file in os.listdir():
        base,ext = os.path.splitext(file)
        if ext == '.psf':
            #print(base+ext)
            fpsf = open(base+ext,'r')
            lines = fpsf.readlines()
            fpsf.close()
            #print(lines[3].split()[5])
            tmp = str(lines[3].split()[5])
            tmp = float(tmp)
            nelc = int(tmp)
            #print(nelc)
            elm = lines[0].split()[0]
            #print(elm)
            dnelc[elm] = nelc
    tot_nelc = 0
    for atom in cell:
        tot_nelc = tot_nelc + dnelc[atom.symbol]
    #print(tot_nelc)
    f = open('./siesta.fdf','r')
    lines = f.readlines()
    f.close()
    f = open('siesta.win','w')
    f.write('num_wann            = '+str(int(tot_nelc/2))+'\n')
    f.write('num_bands           = '+str(int(tot_nelc/2+bandscale))+'\n')
    f.write('translate_home_cell = true'+'\n')
    f.write('guiding_centres     = true'+'\n')
    f.write('write_xyz           = True'+'\n')
    f.write('conv_tol            = 1e-10'+'\n')
    f.write('conv_window         = 5'+'\n')
    f.write('dis_num_iter        = 100000'+'\n')
    f.write('iprint              = 2'+'\n')
    f.write('num_iter            = 10000'+'\n')
    f.write('mp_grid             : 1 1 1'+'\n')
    f.write('gamma_only          = true'+'\n')
    f.write('begin kpoints'+'\n')
    f.write('  0.0000  0.0000  0.0000'+'\n')
    f.write('end kpoints'+'\n')
    f.write('begin atoms_cart'+'\n')
    f.write('Ang'+'\n')
    for atom in cell:
        f.write('  '+atom.symbol+'      '+str(atom.position[0])+'   '+str(atom.position[1])+'   '+str(atom.position[2])+'\n')
    f.write('end atoms_cart'+'\n')
    f.write('begin projections'+'\n')
    f.write('  random'+'\n')
    f.write('end projections'+'\n')
    f.write('begin unit_cell_cart'+'\n')
    f.write('Ang'+'\n')
    wflg = 0
    for line in lines:
        if wflg == 0 and "block LatticeVectors" in line:
            wflg = 1
        elif wflg == 1 and "endblock LatticeVectors" in line:
            wflg = 0
        elif wflg == 1:
            f.write(line)
    #lat = cell.get_cell()
    #f.write('  '+str(lat[0][0])+'     '+str(lat[1][0])+'     '+str(lat[2][0])+'\n')
    #f.write('  '+str(lat[0][1])+'     '+str(lat[1][1])+'     '+str(lat[2][1])+'\n')
    #f.write('  '+str(lat[0][2])+'     '+str(lat[1][2])+'     '+str(lat[2][2])+'\n')
    f.write('end_unit_cell_cart'+'\n')
    f.close()
    
    f = open('siesta.fdf','a')
    f.write('\n')
    f.write('Siesta2Wannier90.NumberOfBands   '+str(int(tot_nelc/2+bandscale))+'\n')
    f.close()

def get_valence(dir_):
    import os
    os.chdir(dir_)
    dnelc = {}
    for file in os.listdir():
        base,ext = os.path.splitext(file)
        if ext == '.psf':
            #print(base+ext)
            fpsf = open(base+ext,'r')
            lines = fpsf.readlines()
            fpsf.close()
            #print(lines[3].split()[5])
            tmp = str(lines[3].split()[5])
            tmp = float(tmp)
            nelc = int(tmp)
            #print(nelc)
            elm = lines[0].split()[0]
            #print(elm)
            dnelc[elm] = nelc
    return dnelc
