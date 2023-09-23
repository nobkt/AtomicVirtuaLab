def mk_packmol_random(mollist,x_box,y_box,z_hox):
    f = open('packmol.inp','w')
    f.write('tolerance 3.0'+'\n')
    f.write('seed -1'+'\n')
    #f.write('randominitialpoint'+'\n')
    f.write('filetype xyz'+'\n')
    f.write('output system.xyz'+'\n')
    f.write('\n')
    for mol in mollist:
        f.write('structure '+str(mol)+'.xyz'+'\n')
        f.write('  number '+str(mollist[mol])+'\n')
        f.write('  inside box 1.0 1.0 1.0 '+str(x_box-1.0)+' '+str(y_box-1.0)+' '+str(z_hox-1.0)+'\n')
        f.write('end structure'+'\n')
        f.write('\n')
    f.close()

def mk_packmol_slab_random(slab,molboxlist,shift=0.0):
    f = open('packmol.inp','w')
    f.write('tolerance 2.0'+'\n')
    f.write('filetype xyz'+'\n')
    f.write('output system.xyz'+'\n')
    f.write('\n')
    f.write('structure '+str(slab)+'.xyz'+'\n')
    f.write(' number 1'+'\n')
    #f.write(' center'+'\n')
    f.write(' fixed 0. 0. 0. 0. 0. 0.'+'\n')
    f.write('end structure'+'\n')
    f.write('\n')
    for mol in molboxlist:
        f.write('structure '+str(molboxlist[mol]['mol'])+'.xyz'+'\n')
        f.write('  number '+str(molboxlist[mol]['num'])+'\n')
        f.write('  inside box '
                +str(molboxlist[mol]['lx'][0]+shift)+' '
                +str(molboxlist[mol]['ly'][0]+shift)+' '
                +str(molboxlist[mol]['lz'][0]+shift)+' '
                +str(molboxlist[mol]['lx'][1]-shift)+' '
                +str(molboxlist[mol]['ly'][1]-shift)+' '
                +str(molboxlist[mol]['lz'][1]-shift)
                +'\n')
        f.write('end structure'+'\n')
        f.write('\n')
    f.close()
        