def qe2dp(dir_,ext):
    from ase.io import read, write
    from ase.io import Trajectory
    import dpdata
    import pathlib
    import os
    #os.system('cp -r '+dir_+' ./tmp')
    #p = pathlib.Path('./tmp')
    p = pathlib.Path(dir_)
    traj = Trajectory('deepmd.traj',mode='a')
    for f in p.glob('**/*.'+str(ext)):
        cell = read(f,index=':')
        for atoms in cell:
            traj.write(atoms=atoms)
    traj.close()
    d_ase=dpdata.MultiSystems.from_file(file_name='./deepmd.traj', fmt='ase/structure')
    d_ase.to_deepmd_raw('deepmd')
    d_ase.to_deepmd_npy('deepmd')
    #os.system('rm deepmd.traj')
    #os.system('rm -rf tmp')

def siesta2dp():
    import dpdata
    d_siesta=dpdata.MultiSystems.from_dir(dir_name='./*',file_name='log_siesta', fmt='siesta/aimd_output')
    print(d_siesta)
    d_siesta.to_deepmd_raw('deepmd')
    d_siesta.to_deepmd_npy('deepmd')
    #os.system('rm deepmd.traj')
    #os.system('rm -rf tmp')

def get_deepmd_list(dir_):
    import os
    dp_list=os.listdir(path=dir_)
    return dp_list

def wt_deepmd_json(dir_,dp_list,rcut,num_steps,decay_steps=5000,prec='low'):
    import re
    import os
    if prec == 'low':
        scale = 1.0
    elif prec == 'medium':
        scale = 1.5
    elif prec == 'high':
        scale = 2.0
    symbols = re.split('[0-9]',dp_list[0])
    symbols = [l for l in symbols if l != '']
    command = 'dp neighbor-stat -s '+dir_+' -r '+str(rcut)+' -t'
    for symbol in symbols:
        command = command + ' '+symbol
    command = command + ' > neigh-stat.log 2>&1'
    os.system(command)
    f = open('neigh-stat.log','r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if 'DEEPMD INFO    max_nbor_size:' in line:
            tmp = re.findall(r"\d+",line)
            max_nbor_size = [float(n)*scale for n in tmp]
    systems=[]
    for l in dp_list:
        systems.append(str(dir_+'/'+str(l)+'/'))
    f = open('./train.json','w')
    f.write('{'+'\n')
    f.write('    "_comment": " model parameters",'+'\n')
    f.write('    "model": {'+'\n')
    f.write('    "type_map":     ')
    lstr = str(symbols)
    f.write(str(lstr.replace("'",'"')))
    f.write(','+'\n')
    f.write('    "descriptor" :{'+'\n')
    f.write('        "type":             "se_e2_a",'+'\n')
    f.write('        "sel":              [')
    nlen = len(max_nbor_size)
    i = 1
    for size in max_nbor_size:
        if i == nlen:
            f.write(' '+str(int(size)))
        else:
            f.write(' '+str(int(size))+',')
        i = i + 1
    f.write('],'+'\n')
    f.write('        "rcut_smth":        0.50,'+'\n')
    f.write('        "rcut":             ')
    f.write(str(rcut))
    f.write(','+'\n')
    f.write('        "neuron":           [25, 50, 100],'+'\n')
    f.write('        "resnet_dt":        false,'+'\n')
    f.write('        "axis_neuron":      16,'+'\n')
    f.write('        "seed":             1,'+'\n')
    f.write('        "_comment":         " thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "fitting_net" : {'+'\n')
    f.write('        "neuron":           [240, 240, 240],'+'\n')
    f.write('        "resnet_dt":        true,'+'\n')
    f.write('        "seed":             1,'+'\n')
    f.write('        "_comment":         " thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "_comment":     " thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "learning_rate" :{'+'\n')
    f.write('    "type":         "exp",'+'\n')
    f.write('    "decay_steps":  '+str(decay_steps)+','+'\n')
    f.write('    "start_lr":     0.001,  '+'\n')
    f.write('    "stop_lr":      3.51e-8,'+'\n')
    f.write('    "_comment":     "thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "loss" :{'+'\n')
    f.write('    "type":         "ener",'+'\n')
    f.write('    "start_pref_e": 0.02,'+'\n')
    f.write('    "limit_pref_e": 1,'+'\n')
    f.write('    "start_pref_f": 1000,'+'\n')
    f.write('    "limit_pref_f": 1,'+'\n')
    f.write('    "start_pref_v": 0.02,'+'\n')
    f.write('    "limit_pref_v": 1,'+'\n')
    f.write('    "_comment":     " thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "training" : {'+'\n')
    f.write('    "training_data": {'+'\n')
    f.write('        "systems":          ')
    lstr = str(systems)
    f.write(str(lstr.replace("'",'"')))
    f.write(','+'\n')
    f.write('        "batch_size":       "auto",'+'\n')
    f.write('        "_comment":         "thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "validation_data":{'+'\n')
    f.write('        "systems":          ')
    lstr = str(systems)
    f.write(str(lstr.replace("'",'"')))
    f.write(','+'\n')
    f.write('        "batch_size":       1,'+'\n')
    f.write('        "numb_btch":        3,'+'\n')
    f.write('        "_comment":         "thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "numb_steps":   ')
    f.write(str(num_steps))
    f.write(','+'\n')
    f.write('    "seed":         10,'+'\n')
    f.write('    "disp_file":    "lcurve.out",'+'\n')
    f.write('    "disp_freq":    100,'+'\n')
    f.write('    "save_freq":    1000,'+'\n')
    f.write('    "_comment":     "thats all"'+'\n')
    f.write('    },'+'\n')
    f.write('    "_comment":             "thats all"'+'\n')
    f.write('}'+'\n')
    f.close()



