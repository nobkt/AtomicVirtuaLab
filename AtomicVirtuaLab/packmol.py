def mk_packmol_random(mollist,x_box,y_box,z_hox):
    f = open('packmol.inp','w')
    f.write('tolerance 2.0'+'\n')
    f.write('filetype xyz'+'\n')
    f.write('output system.xyz'+'\n')
    f.write('\n')
    for mol in mollist:
        f.write('structure '+str(mol)+'.xyz'+'\n')
        f.write('  number '+str(mollist[mol])+'\n')
        f.write('  inside box 0.0 0.0 0.0 '+str(x_box)+' '+str(y_box)+' '+str(z_hox)+'\n')
        f.write(' end structure'+'\n')
        f.write('\n')
    f.close()
        