def molinit(mols):
    d_mols={}
    for mol in mols:
        for mid in range(mols[mol]):
            d_mols[mid+1] = mol
    return d_mols

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

def distance_list(fdata,symbols,a,mol1,b,mol2,d_reaction,sigma):
    from ase.io import read
    from ase.data import atomic_numbers
    Z_of_type={}
    itype=1
    for symbol in symbols:
        Z_of_type[itype] = atomic_numbers[symbol]
        itype = itype + 1
    cell = read(fdata,format='lammps-data',style='full',sort_by_id=True,Z_of_type=Z_of_type)
    
    alist=[]
    for mol in d_reaction:
        if mol1 == mol or mol1+'_' in mol:
            for aid in d_reaction[mol][a]:
                alist.append(aid)
    blist=[]
    for mol in d_reaction:
        if mol2 == mol or mol2+'_' in mol:
            for bid in d_reaction[mol][b]:
                blist.append(bid)
    dis_list=[]
    for aid in alist:
        for bid in blist:
            if aid[1] != bid[1]:
                dis = cell.get_distance(aid[0],bid[0],mic=True,vector=False)
                if dis < sigma:
                    dis_list.append([aid[0],bid[0],dis])
    return dis_list
    



            

                


        
    