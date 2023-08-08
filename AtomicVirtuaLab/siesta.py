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
                             'MD.Target-pressure':str(press)+' GPa',
                             'MD.Target-temperature':str(temp)+' K',
                             'MD.Initial-temperature':str(temp)+' K',
                             'MD.Initial.Time.Step':1,
                             'MD.Final.Time.Step':int(nstep),
                             'MD.Length.Time.Step':str(dt)+' fs',
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':False,
                             'WriteMullikenPop':1,
                             'WriteHirshfeldPop':True,
                             'WriteVoronoiPop':True,
                             'PartialChargesAtEveryGeometry':True,
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
                             'SCFMustConverge':False,
                             'WriteMullikenPop':1,
                             'WriteHirshfeldPop':True,
                             'WriteVoronoiPop':True,
                             'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_optimize(cell,xc,basis_set,mesh_cutoff,kpts,nstep,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized'):
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
                             'MD.NumCGsteps':int(nstep),
                             'SolutionMethod':str(SolutionMethod),
                             'SCFMustConverge':False,
                             'WriteMullikenPop':1,
                             'WriteHirshfeldPop':True,
                             'WriteVoronoiPop':True,
                             'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')

def mk_siesta_input_scf(cell,xc,basis_set,mesh_cutoff,kpts,pseudo_path,SolutionMethod='diagon',MaxSCFIterations=2000,spin='non-polarized'):
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
                             'SCFMustConverge':False,
                             'WriteMullikenPop':1,
                             'WriteHirshfeldPop':True,
                             'WriteVoronoiPop':True,
                             'PartialChargesAtEveryGeometry':True,
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
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
                             'SCFMustConverge':False,
                             'WriteMullikenPop':1,
                             'WriteHirshfeldPop':True,
                             'WriteVoronoiPop':True,
                             'PartialChargesAtEveryGeometry':True,
                             'ExternalElectricField':[str(ex)+' '+str(ey)+' '+str(ez)+' V/Ang'],
                             'SCF.Mixers':['broyden'],
                             'SCF.Mixer.broyden':[('method','broyden'),('weight',0.01),('weight.linear',0.005)]
                             }
              ).write_input(cell,'siesta')


