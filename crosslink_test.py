from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.moltemplate import cp_rdlt,mklt,mk_system_lt,get_chemical_symbols,mklt_fr_mol,rd_lt,ex_lt
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_npt_compress_input_fr_moltemplate
from AtomicVirtuaLab.reaction import set_reac_type, molinit, distance_list
from rdkit import Chem
from ase.io import read
from ase.visualize import view
import os
import shutil
import sys

os.makedirs('./crosslink_test',exist_ok=True)
os.chdir('./crosslink_test')


# 計算セルの長さ(Å)
xbox=100.0
ybox=100.0
zbox=100.0

# モノマーの名前と計算セル中のモノマー分子の数
mols={
    'mol1':16
    }

d_mols = molinit(mols)

# モノマーの名前とSMILES式のリスト
smiles={
    'mol1':"OC(COC(=O)C=C)COC(=O)C=C"
    }

reactions={
    'reac1':[
        ['$atom:C6','@atom:82'],
        ['$atom:C7','@atom:81L'],
        ['$atom:H18','@atom:85'],
        ['$atom:H19','@atom:85LCH2'],
        ['$atom:H20','@atom:85LCH2']
    ],
    'reac2':[
        ['$atom:C12','@atom:82'],
        ['$atom:C13','@atom:81L'],
        ['$atom:H23','@atom:85'],
        ['$atom:H24','@atom:85LCH2'],
        ['$atom:H25','@atom:85LCH2']
    ]
}

reactionids = {
    'mol1':[
        [7,'reac1'],
        [8,'reac1'],
        [13,'reac2'],
        [14,'reac2']
    ]
}

# モノマー分子のsmiles式
os.makedirs('rdltfiles',exist_ok=True)
os.chdir('rdltfiles')
rdltpath = os.getcwd()
mollist={}
bondlist={}
for mol in smiles:
    smi = smiles[mol]
    molname = mol
    # smiles式を読み込む(モノマー)
    smiles2xyz(smi,molname,True,smarts=False,userandom=False)
    monomer = read(molname+'.mol')
    view(monomer)
    cp_rdlt()
    mklt(smi,molname,random=False)
    mollist0,bondlist0 = rd_lt(rdltpath,molname)
    for mol0 in mollist0:
        mollist[mol0] = mollist0[mol0].copy()
    for mol0 in bondlist0:
        bondlist[mol0] = bondlist0[mol0].copy()
os.chdir('../')


# 重合分子のsmiles式
os.makedirs('rdltfiles_poly_test',exist_ok=True)
os.chdir('rdltfiles_poly_test')
smi_poly = "CCC(O)OCC(O)COC(O)CCCC(CCC(=O)OCC(O)COC(=O)CC)C(=O)OCC(O)COC(=O)C(CCCC(=O)OCC(O)COC(=O)CC)CCC(=O)OCC(O)COC(=O)CC"
# smiles式を読み込む(重合)
smiles2xyz(smi_poly,molname+'_poly_test',True,smarts=False,userandom=True)
poly = read(molname+'_poly_test.xyz')
view(poly)
cp_rdlt()
mklt(smi_poly,molname+'_poly_test',random=True)
os.chdir('../')

#mollist = ex_lt(mollist,'mol1',reactions['reac1'])

# 反応後の力場確認

# packmolで分子を充填
os.makedirs('packmolfiles',exist_ok=True)
os.chdir('packmolfiles')
packmolpath = os.getcwd()
for mol in smiles:
    molname = mol
    smi = smiles[mol]
    smiles2xyz(smi,molname,True,smarts=False,userandom=False)
mk_packmol_random(mols,xbox,ybox,zbox)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
os.chdir('..')

# moltemplateの実行(LAMMPS力場ファイル作成)
os.makedirs('moltemplate',exist_ok=True)
os.chdir('moltemplate')
moltemplatepath = os.getcwd()
shutil.copy(packmolpath+'/system.xyz','./')
for mol in smiles:
    molname=mol
    shutil.copy(rdltpath+'/'+str(molname)+'.lt','./')
mk_system_lt(mols,xbox,ybox,zbox)
#os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
#os.system('cleanup_moltemplate.sh')
symbols = get_chemical_symbols('system.data')
d_reaction = set_reac_type('./system.data',symbols,reactionids,d_mols)
os.chdir('../')

sigma = 3.55
dis_list = distance_list(moltemplatepath+'/system.data',symbols,7,'mol1',7,'mol1',d_reaction,sigma)


"""
# NPTで圧縮
os.makedirs('npt_compress',exist_ok=True)
os.chdir('npt_compress')
shutil.copy(moltemplatepath+'/system.data','./')
shutil.copy(moltemplatepath+'/system.in.charges','./')
shutil.copy(moltemplatepath+'/system.in.init','./')
shutil.copy(moltemplatepath+'/system.in.settings','./')
#### TO DO：NPTの圧縮条件を初期条件へ ####
mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,2000,2000,10,10000,400,100,10000,400,200000,300,1.0,'iso',12345,False,qeq=True)
os.system('mpirun -np 4 lmp -in input_npt.lmp 1> log_lmp 2> err_lmp')
os.chdir('../')
"""



