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

def mk_system_lt(d_mols,d_bonds,x_box,y_box,z_box):
    f = open('system.lt','w')
    mollist=[]
    for imol in d_mols:
        if d_mols[imol] in mollist:
            continue
        else:
            mollist.append(d_mols[imol])
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
    for imol in d_mols:
        f.write('mol'+str(imol)+' = new '+d_mols[imol]+' ['+str(1)+']'+'\n')
    f.write('\n')
    if len(d_bonds) != 0:
        f.write("write('Data Bond List') {"+'\n')
        for nbond in d_bonds:
            f.write('  $bond:b'+str(nbond)+' '+d_bonds[nbond][0]+' '+d_bonds[nbond][1]+'\n')
        f.write('  }'+'\n')
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
        elif 'Data Bond List' in line and flg_bonds == False:
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

def ex_lt(path,molname,mollist):
    f = open(path+'/'+molname+'.lt','r')
    lines = f.readlines()
    f.close()
    flg_atoms = False
    ii = 0
    for line in lines:
        if 'Data Atoms' in line and flg_atoms == False:
            flg_atoms = True
        elif 'inherits' in line:
            line0 = line.split()
            if line0[0] == molname:
                for newmol in mollist:
                    lines[ii] = lines[ii].replace(line0[0],newmol)
        elif flg_atoms == True and '}' in line:
            flg_atoms = False
        elif flg_atoms == True:
            line0 = line.split()
            for newmol in mollist:
                for atype in mollist[newmol]:
                    if atype[0] == line0[0]:
                        lines[ii] = lines[ii].replace(line0[2],atype[1])
        ii = ii + 1
    for newmol in mollist:
        f = open(path+'/'+newmol+'.lt','w')
        for line in lines:
            f.write(line)
    f.close()

  


            
            
            
    
    
    
                