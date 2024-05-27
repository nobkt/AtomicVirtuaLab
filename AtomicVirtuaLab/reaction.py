from AtomicVirtuaLab.moltemplate import rd_lt,ex_lt
import AtomicVirtuaLab.globalv as g

def molinit(mols):
    d_mols={}
    for mol in mols:
        for mid in range(mols[mol]):
            d_mols[mid+1] = mol
    return d_mols

def set_reac_atom(fdata,symbols,reactionids,d_mols):
    from ase.io import read
    from ase.data import atomic_numbers
    Z_of_type={}
    itype=1
    for symbol in symbols:
        Z_of_type[itype] = atomic_numbers[symbol]
        itype = itype + 1
    cell = read(fdata,format='lammps-data',style='full',sort_by_id=True,Z_of_type=Z_of_type)
    d_reacatoms={}
    for mol in reactionids:
        for rid in reactionids[mol]:
            natom = 0
            for imol in d_mols:
                if mol == d_mols[imol] or mol+'_' in d_mols[imol]:
                    d_reacatoms[rid[0]+natom]={}
                    d_reacatoms[rid[0]+natom]['mol-id']=imol
                    d_reacatoms[rid[0]+natom]['molname']=d_mols[imol]
                    d_reacatoms[rid[0]+natom]['reacmol']=mol
                    d_reacatoms[rid[0]+natom]['reacid']=rid[0]
                    d_reacatoms[rid[0]+natom]['reactype']=rid[1]
                    d_reacatoms[rid[0]+natom]['reacatom']=rid[2]
                natom = natom + cell.arrays['mol-id'].tolist().count(imol)
    return d_reacatoms

def set_reac_type(fdata,symbols,reactionids,d_mols):
    from ase.io import read
    from ase.data import atomic_numbers
    Z_of_type={}
    itype=1
    for symbol in symbols:
        Z_of_type[itype] = atomic_numbers[symbol]
        itype = itype + 1
    cell = read(fdata,format='lammps-data',style='full',sort_by_id=True,Z_of_type=Z_of_type)
    #cell.arrays=['id']
    #cell.arrays['type']
    #cell.arrays['mol-id']
    d_reaction={}
    for mol in reactionids:
        d_reaction[mol]={}
        for rid in reactionids[mol]:
            d_reaction[mol][rid[0]]=[]
            natom = 0
            for imol in d_mols:
                if mol == d_mols[imol]:
                    d_reaction[mol][rid[0]].append([rid[0]+natom,imol,rid[1]])
                natom = natom + cell.arrays['mol-id'].tolist().count(imol)
    #for mol in d_reaction:
    #    print('mol: ',mol)
    #    for rid in d_reaction[mol]:
    #        print('rid:' ,rid)
    #        for atomid in d_reaction[mol][rid]:
    #            print(atomid)
    return d_reaction

def distance_list(fdata,symbols,d_reacatoms,sigma):
    from ase.io import read
    from ase.data import atomic_numbers
    Z_of_type={}
    itype=1
    for symbol in symbols:
        Z_of_type[itype] = atomic_numbers[symbol]
        itype = itype + 1
    cell = read(fdata,format='lammps-data',style='full',sort_by_id=True,Z_of_type=Z_of_type)
    #
    i = 0
    dis_list=[]
    for ia in d_reacatoms:
        j = 0
        for ib in d_reacatoms:
            if j <= i:
                j = j + 1
                continue
            elif d_reacatoms[ia]['mol-id'] == d_reacatoms[ib]['mol-id']:
                j = j + 1
                continue
            else:
                j = j + 1
                d = cell.get_distance(ia,ib,mic=True,vector=False)
                if d <= sigma:
                    dis_list.append([ia,ib,d])                
        i = i + 1
    return dis_list
    

def reac_test(a,d_reacatoms,d_mols,reactions):
    molname = d_mols[d_reacatoms[a]['mol-id']]
    # Change molname in d_mols
    #if not '_'+d_reacatoms[a]['reactype'] in d_mols[d_reacatoms[a]['mol-id']]:
    #    d_mols[d_reacatoms[a]['mol-id']] = d_mols[d_reacatoms[a]['mol-id']]+'_'+d_reacatoms[a]['reactype']
    d_mols[d_reacatoms[a]['mol-id']] = d_mols[d_reacatoms[a]['mol-id']]+'_'+d_reacatoms[a]['reactype']
    # Change force field
    mollist0,bondlist0 = rd_lt(g.rdltpath,molname)
    mollist0[d_mols[d_reacatoms[a]['mol-id']]] = mollist0[molname]
    del mollist0[molname]
    for mol in mollist0:
        ii = 0
        for atype in mollist0[mol]:
            for rtype in reactions[d_reacatoms[a]['reactype']]:
                #print(rtype)
                if rtype[0] == atype[0]:
                    mollist0[mol][ii][1] = rtype[1]
            ii = ii + 1
    ex_lt(g.rdltpath,molname,mollist0)
    return d_mols


            

                


        
    