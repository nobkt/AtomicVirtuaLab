def cp_rdlt():
    import shutil
    import os
    # rdltをコピーする
    shutil.copy(os.environ.get('RDLTPATH')+'/rdlt.py', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/STaGE_opls_tomoltemplate_opls.txt', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/rdlt.py', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/opls_lt_dict.pkl', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/opls_lt.fdefn', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/lopls_tomoltemplate.txt', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/lopls_lt_dict.pkl', './')
    shutil.copy(os.environ.get('RDLTPATH')+'/lopls_lt.fdefn', './')

def mklt(smiles,molname,random=False):
    import os
    if random:
        os.system('python3 rdlt.py --smi "'+str(smiles)+'" --random True -n '+str(molname)+' -l > '+str(molname)+'.lt')
    else:
        os.system('python3 rdlt.py --smi "'+str(smiles)+'" -n '+str(molname)+' -l > '+str(molname)+'.lt')

def mklt_fr_mol(molfile,molname):
    import os
    os.system('python3 rdlt.py --smi "'+str(molfile)+'" -n '+str(molname)+' > '+str(molname)+'.lt')

def mk_system_lt(mollist,x_box,y_box,z_box):
    f = open('system.lt','w')
    for mol in mollist:
        f.write('import "'+mol+'.lt"'+'\n')
    f.write('\n')
    f.write('write_once("Data Boundary") {'+'\n')
    f.write('  0.0 '+str(x_box)+' xlo xhi'+'\n')
    f.write('  0.0 '+str(y_box)+' ylo yhi'+'\n')
    f.write('  0.0 '+str(z_box)+' zlo zhi'+'\n')
    f.write('}'+'\n')
    f.write('\n')
    for mol in mollist:
        f.write(mol+' = new '+mol+' ['+str(mollist[mol])+']'+'\n')
    f.close()

def get_chemical_symbols(lmpdata):
    from ase.data import chemical_symbols, atomic_masses
    f = open(lmpdata,'r')
    lines = f.readlines()
    f.close()
    rdflg=0
    count=0
    masses=[]
    for line in lines:
        line = line.split()
        if rdflg == 0 and len(line) != 0 and line[0] == 'Masses':
            rdflg = 1
        elif rdflg == 1 and len(line) == 0 and count < 1:
            count = count + 1
        elif rdflg == 1 and len(line) == 0 and count == 1:
            rdflg = 0
        elif rdflg == 1:
            masses.append(float(line[1]))
    symbols=[]
    for mass in masses:
        tmp=[]
        for m in list(atomic_masses):
            tmp.append(abs(float(m)-mass))
        symbols.append(chemical_symbols[tmp.index(min(tmp))])
    return symbols
    
    
            
            
            
    
    
    
                