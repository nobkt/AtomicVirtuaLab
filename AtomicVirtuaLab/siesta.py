def mk_siesta_input_npt(cell,xc,basis_set,mesh_cutoff,kpts,temp,press,nstep,pseudo_path,SolutionMethod='diagon',MixingWeight=0.1,MaxSCFIterations=2000,dt=0.5,spin='non-polarized'):
    from ase.calculators.siesta import Siesta
    from ase.units import Ry
    
    calc = Siesta(label='siesta',
              xc=str(xc),
              mesh_cutoff=float(mesh_cutoff)*Ry,
              basis_set=str(basis_set),
              kpts=kpts,
              spin=spin,
              pseudo_path=pseudo_path,
              fdf_arguments={'DM.MixingWeight':float(MixingWeight),
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
                             'SolutionMethod':str(SolutionMethod)
                             },
              ).write_input(cell,'siesta')

