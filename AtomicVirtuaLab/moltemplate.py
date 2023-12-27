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

def mklt_fr_mol(molfile,molname,random=False):
    import os
    if random:
        os.system("python3 rdlt.py --molfile "+str(molfile)+" --random True -n '+str(molname)+' -l > '+str(molname)+'.lt")
    else:
        os.system("python3 rdlt.py --molfile "+str(molfile)+" -n '+str(molname)+' -l > '+str(molname)+'.lt")

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
    #for mol in mollist:
    #    f.write(mol+' = new '+mol+' ['+str(mollist[mol])+']'+'\n')
    for mol in mollist:
        for imol in range(mollist[mol]):
            f.write('mol'+str(imol+1)+' = new '+mol+' ['+str(1)+']'+'\n')
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

def rd_lt(path,molname):
    import sys
    f = open(path+'/'+molname+'.lt','r')
    lines = f.readlines()
    f.close()
    flg_atoms = False
    flg_bonds = False
    mollist={}
    mollist[molname] = []
    bondlist={}
    bondlist[molname] = []
    for line in lines:
        if 'Data Atoms' in line and flg_atoms == False:
            flg_atoms = True
        if 'Data Bond List' in line and flg_bonds == False:
            flg_bonds = True
        elif flg_atoms == True and '}' in line:
            flg_atoms = False
        elif flg_bonds == True and '}' in line:
            flg_bonds = False
        elif flg_atoms == True:
            line = line.split()
            mollist[molname].append([line[0],line[2]])
        elif flg_bonds == True:
            line = line.split()
            bondlist[molname].append([line[0],line[1],line[2]])
    return mollist, bondlist

def ex_lt(mollist,molname,reaction):
    rid = 1
    if 'r'+str(rid)+'_' in molname:
        rid = rid + 1
    mol_new = 'r'+str(rid)+'_'+molname
    mollist[mol_new] = mollist[molname].copy()
    for l_reac in reaction:
        id=0
        for mols in mollist[mol_new]:
            if mols[0] == l_reac[0]:
                mollist[mol_new][id][1] = l_reac[1]
            id = id + 1
    return mollist
                
        

            
            
    
    
    
            
            
            
    
    
    
                