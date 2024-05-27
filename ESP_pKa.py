from AtomicVirtuaLab.io import smiles2xyz
from ase.io import read
from ase.io.cube import read_cube,read_cube_data
from ase.visualize import view
from ase.calculators.gaussian import Gaussian
from rdkit.Chem import Draw
from rdkit import Chem
import os
import sys

os.makedirs('./ESP_pKa',exist_ok=True)
os.chdir('./ESP_pKa')


# cubeファイル読み込み
l_xc=['B3LYP','wB97XD','HF']
l_isodens=[0.00002]
l_mesh=[0]
diff=5e-5

"""
os.makedirs('data1',exist_ok=True)
os.chdir('data1')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    for mesh in l_mesh:
        os.makedirs(str(mesh),exist_ok=True)
        os.chdir(str(mesh))
        for isodens in l_isodens:
            os.makedirs(str(isodens),exist_ok=True)
            os.chdir(str(isodens))
            f=open('MESP.csv','w')
            for imol in range(45):
                premax=[]
                premin=[]
                espdata,esp = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data1/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian.cube')
                densdata,dens = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data1/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian_dens.cube')
                for d in zip(espdata,densdata):
                    esp0=[]
                    for dd in zip(d[0],d[1]):
                        for ddd in zip(dd[0],dd[1]):
                            if  ddd[1] > isodens-diff and ddd[1] < isodens+diff:
                                esp0.append(ddd[0])
                            else:
                                continue
                    if len(esp0) != 0:
                        premax.append(max(esp0))
                        premin.append(min(esp0))
                if len(premax) == 0 or len(premin) == 0:
                    continue
                else:
                    f.write(str(imol+1)+','+str(max(premax))+' '+str(min(premin))+'\n')
                    print(xc,isodens,imol+1,max(premax),min(premin))
            f.close()
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')


os.makedirs('data2',exist_ok=True)
os.chdir('data2')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    for mesh in l_mesh:
        os.makedirs(str(mesh),exist_ok=True)
        os.chdir(str(mesh))
        for isodens in l_isodens:
            os.makedirs(str(isodens),exist_ok=True)
            os.chdir(str(isodens))
            f=open('MESP.csv','w')
            for imol in range(22):
                premax=[]
                premin=[]
                espdata,esp = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data2/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian.cube')
                densdata,dens = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data2/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian_dens.cube')
                for d in zip(espdata,densdata):
                    esp0=[]
                    for dd in zip(d[0],d[1]):
                        for ddd in zip(dd[0],dd[1]):
                            if  ddd[1] > isodens-diff and ddd[1] < isodens+diff:
                                esp0.append(ddd[0])
                            else:
                                continue
                    if len(esp0) != 0:
                        premax.append(max(esp0))
                        premin.append(min(esp0))
                if len(premax) == 0 or len(premin) == 0:
                    continue
                else:
                    f.write(str(imol+1)+','+str(max(premax))+' '+str(min(premin))+'\n')
                    print(xc,isodens,imol+1,max(premax),min(premin))
            f.close()
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')


os.makedirs('data3',exist_ok=True)
os.chdir('data3')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    for mesh in l_mesh:
        os.makedirs(str(mesh),exist_ok=True)
        os.chdir(str(mesh))
        for isodens in l_isodens:
            os.makedirs(str(isodens),exist_ok=True)
            os.chdir(str(isodens))
            f=open('MESP.csv','w')
            for imol in range(5):
                premax=[]
                premin=[]
                espdata,esp = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data3/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian.cube')
                densdata,dens = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data3/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian_dens.cube')
                for d in zip(espdata,densdata):
                    esp0=[]
                    for dd in zip(d[0],d[1]):
                        for ddd in zip(dd[0],dd[1]):
                            if  ddd[1] > isodens-diff and ddd[1] < isodens+diff:
                                esp0.append(ddd[0])
                            else:
                                continue
                    if len(esp0) != 0:
                        premax.append(max(esp0))
                        premin.append(min(esp0))
                if len(premax) == 0 or len(premin) == 0:
                    continue
                else:
                    f.write(str(imol+1)+','+str(max(premax))+' '+str(min(premin))+'\n')
                    print(xc,isodens,imol+1,max(premax),min(premin))
            f.close()
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

os.makedirs('data4',exist_ok=True)
os.chdir('data4')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    for mesh in l_mesh:
        os.makedirs(str(mesh),exist_ok=True)
        os.chdir(str(mesh))
        for isodens in l_isodens:
            os.makedirs(str(isodens),exist_ok=True)
            os.chdir(str(isodens))
            f=open('MESP.csv','w')
            for imol in range(1):
                premax=[]
                premin=[]
                espdata,esp = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data3/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian.cube')
                densdata,dens = read_cube_data('/home/A23321P/work/myGaussian/ESP_pKa/data3/'+str(xc)+'/mol'+str(imol+1)+'/'+str(mesh)+'/gaussian_dens.cube')
                for d in zip(espdata,densdata):
                    esp0=[]
                    for dd in zip(d[0],d[1]):
                        for ddd in zip(dd[0],dd[1]):
                            if  ddd[1] > isodens-diff and ddd[1] < isodens+diff:
                                esp0.append(ddd[0])
                            else:
                                continue
                    if len(esp0) != 0:
                        premax.append(max(esp0))
                        premin.append(min(esp0))
                if len(premax) == 0 or len(premin) == 0:
                    continue
                else:
                    f.write(str(imol+1)+','+str(max(premax))+' '+str(min(premin))+'\n')
                    print(xc,isodens,imol+1,max(premax),min(premin))
            f.close()
            os.chdir('../')
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
#
"""


# SMILESから静電ポテンシャル計算
l_xc =[
    "wB97XD",
    "B3LYP",
    "HF"
]

l_smiles1 = [
    "[O-][N+](=C)C1=CC(S)=CC=C1",
    "SC1=CNC(=O)NC1=O",
    "SC1=CC=C(C=O)C=C1",
    "SC1=CC(Cl)=CC=C1",
    "SC1=CC=C(Cl)C=C1",
    "SC1=CC=CC=C1",
    "CC1=C(S)C=CC=C1",
    "CC1=CC=CC(S)=C1",
    "COC1=CC=C(S)C=C1",
    "CC1=CC=C(S)C=C1",
    "FC(F)(F)CS",
    "CC(S)=C",
    "CCOC(=O)CS",
    "COC(=O)CS",
    "CC(C)C(=O)NCC(O)=O",
    "OCC(S)CS",
    "CCOCCS",
    "SCC1=CC=CC=C1",
    "OCC(O)CS",
    "OCCS",
    "CC(C)(S)CO",
    "SCC=C",
    "CS",
    "CCS",
    "CCCCS",
    "CC(C)S",
    "CC(O)(S)CS",
    "N[C@@H](CS)C(O)=O",
    "OC(=O)CS",
    "NCCS",
    "OCC(S)CS",
    "CC(C)(C)S",
    "CCC(C)(C)S",
    "SCCS",
    "CC(S)C([O-])=O",
    "SCC1=NC=CC=C1",
    "SC1CCCCC1",
    "COC1=CC=CC(S)=C1",
    "SC1=CC=C(Br)C=C1",
    "[O-][N+](=O)C1=CC=C(S)C=C1",
    "OC(=O)CCS",
    "OCCCS",
    "COC(=O)CCS",
    "SCCCCS",
    "SCCNC(=O)NCCS"
]

l_smiles2 = [
    "CC(O)=O",
    "CCC(=O)O",
    "CCCC(=O)O",
    "CCCCC(=O)O",
    "CCCCCC(O)=O",
    "C(C(=O)O)Cl",
    "C(C(=O)O)Br",
    "C(=O)(C(Cl)(Cl)Cl)O",
    "C1=CC=C(C(=C1)C(=O)O)Cl",
    "CC(CC(=O)O)Cl",
    "C(CC(=O)O)CCl",
    "C=CCC(=O)O",
    "CC(C)C(=O)O",
    "CC(C)(C)C(=O)O",
    "CC(C)CC(=O)O",
    "CCC(C)C(=O)O",
    "CC#CC(=O)O",
    "CC(C(=O)O)Cl",
    "C(CBr)C(=O)O",
    "C(CCl)C(=O)O",
    "C\C=C\C(O)=O",
    "C(=O)O"
]

l_smiles3=[
    "CN(C)C1=CC=NC=C1",
    "OC([O-])=O",
    "[OH-]",
    "CC([O-])=O",
    "CCCCC(CC)C([O-])=O"
]

l_smiles4=[
    "[H][N+](C)(C)C1=CC=NC=C1",
    "[H][N+]1=CC=C(C=C1)N(C)C",
    "[H]OC(O)=O",
    "O",
    "[H]OC(C)=O",
    "[H]OC(=O)C(CC)CCCC"
]

os.makedirs('data1',exist_ok=True)
os.chdir('data1')

for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    imol=1
    for smiles in l_smiles1:
        os.makedirs('mol'+str(imol),exist_ok=True)
        os.chdir('mol'+str(imol))
        Draw.MolToFile(Chem.MolFromSmiles(smiles),'mol'+str(imol)+'.png')
        smiles2xyz(smiles,'mol'+str(imol),True,smarts=False)
        #smiles2xyz(smiles,'test',True,smarts=False,userandom=rnd)
        mol_ = read('mol'+str(imol)+'.xyz')
        #view(mol_)
        if xc == 'PM6':
            Gaussian(label='gaussian',
                nprocshared=16,
                chk='gaussian.chk',
                xc=xc,
                #basis='6-311++G(d,p)',
                scf='qc,tight,maxcycle=9999',
                opt='calcfc,maxcycles=9999',
                freq=''
                ).write_input(mol_)
        else:
            Gaussian(label='gaussian',
                nprocshared=16,
                chk='gaussian.chk',
                xc=xc,
                basis='6-311++G(d,p)',
                scf='qc,tight,maxcycle=9999',
                opt='calcfc,maxcycles=9999',
                freq=''
                ).write_input(mol_)
        imol = imol + 1
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

os.makedirs('data2',exist_ok=True)
os.chdir('data2')


for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    imol=1
    for smiles in l_smiles2:
        os.makedirs('mol'+str(imol),exist_ok=True)
        os.chdir('mol'+str(imol))
        Draw.MolToFile(Chem.MolFromSmiles(smiles),'mol'+str(imol)+'.png')
        smiles2xyz(smiles,'mol'+str(imol),True,smarts=False)
        #smiles2xyz(smiles,'test',True,smarts=False,userandom=rnd)
        mol_ = read('mol'+str(imol)+'.xyz')
        #view(mol_)
        Gaussian(label='gaussian',
            nprocshared=16,
            chk='gaussian.chk',
            xc=xc,
            basis='6-311++G(d,p)',
            scf='qc,tight,maxcycle=9999',
            opt='calcfc,maxcycles=9999',
            freq=''
            ).write_input(mol_)
        imol = imol + 1
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

os.makedirs('data3',exist_ok=True)
os.chdir('data3')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    imol=1
    for smiles in l_smiles3:
        os.makedirs('mol'+str(imol),exist_ok=True)
        os.chdir('mol'+str(imol))
        Draw.MolToFile(Chem.MolFromSmiles(smiles),'mol'+str(imol)+'.png')
        smiles2xyz(smiles,'mol'+str(imol),True,smarts=False)
        #smiles2xyz(smiles,'test',True,smarts=False,userandom=rnd)
        mol_ = read('mol'+str(imol)+'.xyz')
        #view(mol_)
        Gaussian(label='gaussian',
            nprocshared=16,
            chk='gaussian.chk',
            xc=xc,
            basis='6-311++G(d,p)',
            scf='qc,tight,maxcycle=9999',
            opt='calcfc,maxcycles=9999',
            freq=''
            ).write_input(mol_)
        imol = imol + 1
        os.chdir('../')
    os.chdir('../')
os.chdir('../')

os.makedirs('data4',exist_ok=True)
os.chdir('data4')
for xc in l_xc:
    os.makedirs(xc,exist_ok=True)
    os.chdir(xc)
    imol=1
    for smiles in l_smiles4:
        os.makedirs('mol'+str(imol),exist_ok=True)
        os.chdir('mol'+str(imol))
        Draw.MolToFile(Chem.MolFromSmiles(smiles),'mol'+str(imol)+'.png')
        smiles2xyz(smiles,'mol'+str(imol),True,smarts=False)
        #smiles2xyz(smiles,'test',True,smarts=False,userandom=rnd)
        mol_ = read('mol'+str(imol)+'.xyz')
        #view(mol_)
        Gaussian(label='gaussian',
            nprocshared=16,
            chk='gaussian.chk',
            xc=xc,
            basis='6-311++G(d,p)',
            scf='qc,tight,maxcycle=9999',
            opt='calcfc,maxcycles=9999',
            freq=''
            ).write_input(mol_)
        imol = imol + 1
        os.chdir('../')
    os.chdir('../')
os.chdir('../')
# SMILESから静電ポテンシャル計算
