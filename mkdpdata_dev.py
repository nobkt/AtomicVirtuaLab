from AtomicVirtuaLab.deepmd import qe2dp, get_deepmd_list, wt_deepmd_json
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.espresso import mk_qe_input_npt, mk_qe_input_scf, mk_qe_input_nvt
from AtomicVirtuaLab.lammps import mk_npt_input_dpdata, mk_nvt_input_dpdata, mk_npt_slab_input_dpdata
from AtomicVirtuaLab.io import rd_cif, cell2atomlist
from ase.io import read
from ase.visualize import view
from ase.build import make_supercell
from ase.data import atomic_numbers
from ase.units import Rydberg
import os
import sys
import random
import math

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifdir='/home/A23321P/work/myPython/AtomicVirtuaLab/cifs/dpdatas'

ncore=32
ngpu=1

usepbs = True
host='usagi3'
mem=128

dpexedir='/home/modules/applications/gpu/deepmd-kit-2.2.7-cuda10.2/bin'
qeexedir='/home/modules/applications/cpu/qe/7.3'

d_cifs = {
    'MgAl2O4_primitive' : [2,2,2],
    'MgAl2O4_slab100_1' : [1,1,1],
    'MgAl2O4_slab100_2' : [1,1,1],
    'MgAl2O4_slab110_1' : [1,1,1],
    'MgAl2O4_slab110_2' : [1,1,1],
    'MgAl2O4_slab111_1' : [1,1,1],
    'MgAl2O4_slab111_2' : [1,1,1],
    'MgAl2O4_slab111_3' : [1,1,1],
    'MgAl2O4_slab111_4' : [1,1,1]
}


d_input = {
    'MgAl2O4_primitive' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'all'
    },
    'MgAl2O4_slab100_1' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab100_2' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab110_1' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab110_2' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab111_1' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab111_2' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab111_3' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    },
    'MgAl2O4_slab111_4' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [300,1000,2000,3000,4000,5000],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 5000,
        'Tstep' : 100,
        'mode' : 'npt',
        'cell_dofree' : 'c'
    }
}

"""
d_cifs = {
    'Si' : [2,2,2],
    '3C-SiC' : [2,2,2]
}

d_input = {
    'Si' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [100],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 100,
        'Tstep' : 100,
        'mode' : 'npt'
    },
    '3C-SiC' : {
        'ninit' : 500,
        'nout' : 10,
        'kpts' : None,
        'l_temp' : [100],
        'l_press' : [0],
        'Tsta' : 100,
        'Tend' : 100,
        'Tstep' : 100,
        'mode' : 'npt'
    }
}
"""


options={'vdw_corr':'dft-d3','dftd3_version':4}
xc = 'pbesol'

dEaccept=4.336454e-2
f_hi = 0.15
f_lo = 0.05
final_accept = 99.0
maxset = 100

restart = False
ntrain0 = 1
ntraj0 = 151
init_skip = False

#ninit = 500
#nout = 10
#ninit = 10
#nout = 10
#kpts=(2,2,2)
#kpts=None

#l_temp = [300, 1000, 2000, 3000, 4000, 5000]
#l_temp = [300]
#l_temp = [1000]
#l_press = [0]

#Tsta=100
#Tend=3000
#Tsta=100
#Tend=100
#Tstep=100

#os.environ['CUDA_VISIBLE_DEVICES'] = "0"
os.makedirs('mkdpdatas',exist_ok=True)
os.chdir('mkdpdatas')
main_path = os.getcwd()

#
# 1. read CIF
#
d_strc = {}
for fcif in d_cifs:
    cell = rd_cif(g.cifdir+'/'+fcif+'.cif',primitive_cell=False)
    lat = d_cifs[fcif]
    cell = make_supercell(cell,([lat[0],0,0],[0,lat[1],0],[0,0,lat[2]]))
    d_strc[fcif] = cell
#for str in d_str:
#    view(d_str[str])

#
# 2. make initial train datas
#
#trajs=[]
nset=1
d_dppath={}
if restart:
    ftrain = open('train.out','a')
else:
    ftrain = open('train.out','w')
    ftrain.write('ntrain ntry naccept accept_rate'+'\n')
for strc in d_strc:
    os.makedirs(strc,exist_ok=True)
    os.chdir(strc)
    ntraj=1
#
# 2-1. QE NPT
#
    l_temp = d_input[strc]['l_temp']
    l_press = d_input[strc]['l_press']
    ninit = d_input[strc]['ninit']
    kpts = d_input[strc]['kpts']
    mode = d_input[strc]['mode']
    cell_dofree = d_input[strc]['cell_dofree']
    for temp in l_temp:
        for press in l_press:
            os.makedirs('qe',exist_ok=True)
            os.chdir('qe')
            qe_path = os.getcwd()
            if ntraj == 1:
                dppath = os.getcwd()
                d_dppath[strc] = dppath
            if not init_skip:
                os.makedirs('set'+str(nset),exist_ok=True)
                os.chdir('set'+str(nset))
                #if ntraj == 1:
                #    dppath = os.getcwd()
                #    d_dppath[strc] = dppath
                os.makedirs('run',exist_ok=True)
                os.chdir('run')
                cell = d_strc[strc]
                #os.system("echo 'cd $PBS_O_WORKDIR;mpirun -np "+str(ncore)+" pw.x < qe_npt.pwi 1> qe_npt.pwo 2> qe_npt.err' | qsub -l select=host=usagi3:ncpus=64 -o finishjob")
                #while True:
                #    if os.path.isfile('./finishjob'):
                #        break
                #traj = read('qe_npt.pwo',index=':')
                #for atoms in traj:
                #    trajs.append(atoms)
                if mode == 'npt':
                    mk_qe_input_npt(cell,xc,'paw',temp,100,press,dt=0.5,level='SSSP_precision',estep=9999,nstep=ninit,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=kpts,ecut='auto',cell_dofree='all',options=options,cell_dofree=cell_dofree,nspin=False,magElm=None)
                    if usepbs:
                        os.system("echo 'source "+str(dpexedir)+"/activate;source "+str(qeexedir)+"/qe_profile.sh;cd $PBS_O_WORKDIR;mpirun -np "+str(ncore)+" pw.x < qe_npt.pwi 1> qe_npt.pwo 2> qe_npt.err' | qsub -l select=1:mem="+str(mem)+"gb:host="+str(host)+":ncpus="+str(ncore)+" -o finishjob")
                        while True:
                            if os.path.isfile('./finishjob'):
                                os.system('rm finishjob')
                                break
                    else:
                        os.system('mpirun -np '+str(ncore)+' pw.x < qe_npt.pwi 1> qe_npt.pwo 2> qe_npt.err')
                    os.chdir('../')
                    os.system('cp ./run/qe_npt.pwo ./'+str(ntraj)+'.pwo')
                elif mode == 'nvt':
                    #mk_qe_input_nvt(cell,'pbe','paw',temp,100,press,dt=0.5,level='SSSP_precision',estep=9999,nstep=ninit,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=kpts,ecut='auto',cell_dofree='all',options=options,nspin=False,magElm=None)
                    mk_qe_input_nvt(cell,xc,'paw',temp,100,dt=0.5,level='SSSP_precision',estep=9999,nstep=ninit,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=kpts,ecut='auto',options=options,nspin=False)
                    if usepbs:
                        os.system("echo 'source "+str(dpexedir)+"/activate;source "+str(qeexedir)+"/qe_profile.sh;cd $PBS_O_WORKDIR;mpirun -np "+str(ncore)+" pw.x < qe_nvt.pwi 1> qe_nvt.pwo 2> qe_nvt.err' | qsub -l select=1:mem="+str(mem)+"gb:host="+str(host)+":ncpus="+str(ncore)+" -o finishjob")
                        while True:
                            if os.path.isfile('./finishjob'):
                                os.system('rm finishjob')
                                break
                    else:
                        os.system('mpirun -np '+str(ncore)+' pw.x < qe_nvt.pwi 1> qe_nvt.pwo 2> qe_nvt.err')
                    os.chdir('../')
                    os.system('cp ./run/qe_nvt.pwo ./'+str(ntraj)+'.pwo')
                ntraj = ntraj + 1
                os.system('rm -rf run')
                os.chdir('../')
            os.chdir('../')
    os.chdir('../')
    ntraj=1
#
# 2-2. translate to dp format and train
#
os.makedirs('dpdatas',exist_ok=True)
os.chdir('dpdatas')
dp_path = os.getcwd()
if not init_skip:
    os.makedirs('set'+str(nset),exist_ok=True)
    os.chdir('set'+str(nset))
    os.makedirs('run',exist_ok=True)
    os.chdir('run')
    dp_list_path = os.getcwd()
    qe2dp(main_path,'pwo')
    dp_list=get_deepmd_list(dp_list_path+'/deepmd')
    wt_deepmd_json(dp_list_path+'/deepmd',dp_list,6.0,100000,decay_steps=500,prec='high')
    os.system("dp train train.json 1> log_train 2> err_train")
    os.system("dp freeze -o graph.pb 1> log_freeze 2> err_freeze")
    os.chdir('../')
    os.system('cp -r ./run/deepmd .')
    os.system('cp ./run/lcurve.out .')
    os.system('cp ./run/graph.pb .')
    os.system('rm -rf run')
    dp_list_path = os.getcwd()
    os.chdir('../')
os.chdir('../')
os.makedirs('forcefield',exist_ok=True)
os.chdir('forcefield')
forcefield_path = os.getcwd()
if not init_skip:
    os.makedirs('set'+str(nset),exist_ok=True)
    os.chdir('set'+str(nset))
    os.system('cp '+dp_list_path+'/graph.pb .')
    os.system('rm '+dp_list_path+'/graph.pb')
    os.chdir('..')
    init_skip = False
os.chdir('../')
nset = nset + 1  
#
# 3. start train loop
#
#
runflg = 0
for ntrain in range(maxset):
    if runflg == 1:
        nset=1            
        break
    if restart:
        if ntrain <= ntrain0:
            nset = nset + 1
            ntraj = ntraj0
            continue
#
# 3.1 make lammps trajectory
#
    naccpept = 0
    ntry = 0
    for strc in d_strc:
        Tsta = d_input[strc]['Tsta']
        Tend = d_input[strc]['Tend']
        Tstep = d_input[strc]['Tstep']
        nout = d_input[strc]['nout']
        ninit = d_input[strc]['ninit']
        kpts = d_input[strc]['kpts']
        mode = d_input[strc]['mode']
        cell_dofree = d_input[strc]['cell_dofree']
        os.makedirs(strc,exist_ok=True)
        os.chdir(strc)
        strc_path = os.getcwd()
        os.makedirs('lammps',exist_ok=True)
        os.chdir('lammps')
        lammps_path = os.getcwd()
        os.makedirs('set'+str(nset),exist_ok=True)
        os.chdir('set'+str(nset))
        #naccpept = 0
        #ntry = 0
        nlmpset=1
        for temp in range(Tsta,Tend+Tstep,Tstep):
            for press in l_press:
                os.makedirs('lmpset'+str(nlmpset),exist_ok=True)
                os.chdir('lmpset'+str(nlmpset))
                os.makedirs('run',exist_ok=True)
                os.chdir('run')
                cell = d_strc[strc]
                seed=random.randint(10000, 99999)
                #print('seed = ',seed,flush=True)
                if mode == 'npt':
                    if cell_dofree == 'all':
                        mk_npt_input_dpdata(cell,0.00005,100,20000,200,nout,temp,press,seed,mol=False)
                    else:
                        mk_npt_slab_input_dpdata(cell,0.00005,100,20000,200,nout,temp,press,seed,mol=False)
                elif mode == 'nvt':
                    mk_nvt_input_dpdata(cell,0.00005,100,20000,200,nout,temp,seed,mol=False)
                os.system('cp '+forcefield_path+'/set'+str(nset-1)+'/graph.pb .')
                os.system(dpexedir+'/mpirun -np '+str(ngpu)+' lmp -in lammps.lmp 1> log_lmp 2> err_lmp')
                symbols = cell2atomlist(cell)
                ntype=1
                Z_of_type={}
                for symbol in symbols:
                    Z_of_type[ntype] = atomic_numbers[symbol]
                    ntype=ntype+1
                trajs=[]
                for i in range(nout):
                    cell = read('result'+str(i+1)+'.data',format='lammps-data',Z_of_type=Z_of_type,style='atomic')
                    trajs.append(cell)
                f = open('log_lmp','r')
                lines = f.readlines()
                f.close()
                lmp_ene=[]
                for line in lines:
                    if 'Potential Energy:' in line:
                       line = line.split()
                       lmp_ene.append(float(line[4]))
                for i in range(nout):
                    f = open('dump'+str(i+1)+'.lammpstrj.200')
                    lines = f.readlines()
                    f.close()
                    j = 0
                    lmp_force={}
                    for line in lines:
                        if 'ITEM: NUMBER OF ATOMS' in line:
                            line = lines[j+1].split()
                            natom = int(line[0])
                        elif 'ITEM: ATOMS id type fx fy fz' in line:
                            for k in range(natom):
                                line = lines[j+1+k].split()
                                lmp_force[str(line[0])] = [float(line[2]),float(line[3]),float(line[4])]
                        j = j + 1
                os.chdir('../')
                os.system('cp ./run/result* .')
                os.system('cp ./run/log_lmp .')
                os.system('rm -rf run')
                os.chdir('../')
                nlmpset = nlmpset + 1
#
# 3.2. run QE using lammps trajectory's structure
#
                os.chdir(strc_path+'/qe')
                if ntraj == 1:
                    dppath = os.getcwd()
                    d_dppath[strc] = dppath
                os.makedirs('set'+str(nset),exist_ok=True)
                os.chdir('set'+str(nset))
                #if ntraj == 1:
                #    dppath = os.getcwd()
                #    d_dppath[strc] = dppath
                qe_ene=[]
                qe_force={}
                for i in range(nout):
                    os.makedirs('run',exist_ok=True)
                    os.chdir('run')
                    mk_qe_input_scf(trajs[i],xc,'paw',level='SSSP_precision',estep=9999,ecutwfc=25,ecutrho=225,mixing_beta=0.2,kpts=kpts,ecut='auto',tstress=True,nosym=True,options=options,nspin=False)
                    if usepbs:
                        os.system("echo 'source "+str(dpexedir)+"/activate;source "+str(qeexedir)+"/qe_profile.sh;cd $PBS_O_WORKDIR;mpirun -np "+str(ncore)+" pw.x < qe_scf.pwi 1> qe_scf.pwo 2> qe_scf.err' | qsub -l select=1:mem="+str(mem)+"gb:host="+str(host)+":ncpus="+str(ncore)+" -o finishjob")
                        while True:
                            if os.path.isfile('./finishjob'):
                                os.system('rm finishjob')
                                break
                    else:
                        os.system('mpirun -np '+str(ncore)+' pw.x < qe_scf.pwi 1> qe_scf.pwo 2> qe_scf.err')
                    #os.system("echo 'cd $PBS_O_WORKDIR;mpirun -np "+str(ncore)+" pw.x < qe_npt.pwi 1> qe_npt.pwo 2> qe_npt.err' | qsub -l select=host=usagi3:ncpus=64 -o finishjob")
                    #while True:
                    #    if os.path.isfile('./finishjob'):
                    #        break
                    f = open('qe_scf.pwo','r')
                    lines = f.readlines()
                    f.close()
                    for line in lines:
                        if '!    total energy' in line:
                            line = line.split()
                            qe_ene.append(float(line[4]))
                    f = open('qe_scf.pwo','r')
                    lines = f.readlines()
                    f.close()
                    j = 0
                    for line in lines:
                        if 'Forces acting on atoms (cartesian axes, Ry/au):' in line:
                            for k in range(natom):
                                line = lines[j+2+k].split()
                                qe_force[str(line[1])] = [float(line[6])*Rydberg,float(line[7])*Rydberg,float(line[8])*Rydberg]
                        j = j + 1                   
                    #print('E(QE) = ',qe_ene[i]*Rydberg,' ','E(LAMMPS) =',lmp_ene[i])
#
# 3.3. update train data
#
                    sigma = 0.0
                    for nn in range(natom):
                        sigma = sigma + (lmp_force[str(nn+1)][0]-qe_force[str(nn+1)][0])**2
                        sigma = sigma + (lmp_force[str(nn+1)][1]-qe_force[str(nn+1)][1])**2
                        sigma = sigma + (lmp_force[str(nn+1)][2]-qe_force[str(nn+1)][2])**2
                        sigma = sigma/float(3*natom)
                    kappa = 2.0*math.sqrt(sigma)
                    dE = abs(qe_ene[i]*Rydberg-lmp_ene[i])
                    #if dE <= dEaccept:
                    #    os.system('mv qe_scf.pwo ../'+str(ntraj)+'.pwo')
                    #    naccpept = naccpept + 1
                    #    ntraj = ntraj + 1
                    if kappa >= f_lo and kappa <= f_hi:
                        os.system('mv qe_scf.pwo ../'+str(ntraj)+'.pwo')
                        naccpept = naccpept + 1
                        ntraj = ntraj + 1
                    elif kappa < f_lo:
                        naccpept = naccpept + 1
                    ntry = ntry + 1
                    #print(strc,ntrain,ntry,naccpept,dE,flush=True)
                    ftrain.write(strc+' '+str(ntrain)+' '+str(ntry)+' '+str(naccpept)+' '+str(dE)+'\n')
                    ftrain.flush()
                    os.chdir('../')
                    os.system('rm -rf run')
                os.chdir(lammps_path+'/set'+str(nset))
        os.chdir(main_path)
    accept_rate = float(naccpept)/float(ntry)*100.0
    #print(strc,ntrain,ntry,naccpept,accept_rate,flush=True)
    ftrain.write(strc+' '+str(ntrain)+' '+str(ntry)+' '+str(naccpept)+' '+str(accept_rate)+'\n')
    ftrain.flush()
    if accept_rate >= final_accept:
        ftrain.close()
        runflg = 1
#
# 3.3. translate to dp format and train
#
    if runflg == 0:
        os.chdir(dp_path)
        os.makedirs('set'+str(nset),exist_ok=True)
        os.chdir('set'+str(nset))
        os.makedirs('run',exist_ok=True)
        os.chdir('run')
        dp_list_path = os.getcwd()
        qe2dp(main_path,'pwo')
        dp_list=get_deepmd_list(dp_list_path+'/deepmd')
        wt_deepmd_json(dp_list_path+'/deepmd',dp_list,6.0,100000,decay_steps=500,prec='high')
        os.system("dp train train.json 1> log_train 2> err_train")
        os.system("dp freeze -o graph.pb 1> log_freeze 2> err_freeze")
        os.chdir('../')
        os.system('cp -r ./run/deepmd .')
        os.system('cp ./run/lcurve.out .')
        os.system('cp ./run/graph.pb .')
        os.system('rm -rf run')
        dp_list_path = os.getcwd()
        os.chdir('..')
        os.chdir('../')
        os.chdir(forcefield_path)
        os.makedirs('set'+str(nset),exist_ok=True)
        os.chdir('set'+str(nset))
        os.system('cp '+dp_list_path+'/graph.pb .')
        os.system('rm '+dp_list_path+'/graph.pb')
        os.chdir('..')
        os.chdir('../')
        nset = nset + 1
    os.chdir(main_path)
ftrain.close()
#view(trajs)

    



