def mk_nwchem_input_opt(mol,label,xc,basis='6-31G**',chg=0,mult=1):
    from ase.calculators.nwchem import NWChem
    import shutil
    calc = NWChem(
        label=label,
        basispar='spherical',
        basis={
            '*': str(basis)
            },
        charge=chg,
        dft={
            'xc': str(xc),
            'iterations': 9999,
            'semidirect': None,
            'mult': mult
            },
        set={
            'quickguess': True
            },
        driver={'maxiter': 2000,
                'inhess': 1},
        esp={
            'restrain': None
            },
        task='optimize'
        
    ).write_input(mol)
    f = open(label+'.nwi','a')
    print('task esp',file=f)
    f.close()
    f = open(label+'.nwi','r')
    lines = f.readlines()
    f.close()
    f = open(label+'.nwi','w')
    for line in lines:
        if 'permanent_dir' in line:
            continue
        elif 'scratch_dir' in line:
            continue
        else:
            f.write(line)
    f.close()
    shutil.rmtree ('./'+str(label))

    
def mk_nwchem_input_scf(mol,label,xc,basis='6-31G**',chg=0,mult=1):
    from ase.calculators.nwchem import NWChem
    import shutil
    calc = NWChem(
        label=label,
        basispar='spherical',
        basis={
            '*': str(basis)
            },
        charge=chg,
        dft={
            'xc': str(xc),
            'iterations': 9999,
            'semidirect': None,
            'mult': mult
            },
        set={
            'quickguess': True
            },
        driver={'maxiter': 2000
                },
        esp={
            'restrain': None
            },
    ).write_input(mol)
    f = open(label+'.nwi','a')
    print('task esp',file=f)
    f.close()
    f = open(label+'.nwi','r')
    lines = f.readlines()
    f.close()
    f = open(label+'.nwi','w')
    for line in lines:
        if 'permanent_dir' in line:
            continue
        elif 'scratch_dir' in line:
            continue
        else:
            f.write(line)
    f.close()
    shutil.rmtree ('./'+str(label))

    
