"""
from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.moltemplate import cp_rdlt,mklt,mk_system_lt,get_chemical_symbols,mklt_fr_mol,rd_lt,ex_lt
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_npt_compress_input_fr_moltemplate
from AtomicVirtuaLab.reaction import set_reac_type, molinit, distance_list, set_reac_atom, reac_test
import AtomicVirtuaLab.globalv as g
from rdkit import Chem
from ase.io import read
from ase.visualize import view
import os
import shutil
import sys
"""

import os
import random
import sys
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write
from ase.visualize import view
import AtomicVirtuaLab.globalv as g
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_nvt_input_uff_rigid_scale
from AtomicVirtuaLab.build import sortmol


os.makedirs('./crosslink_test',exist_ok=True)
os.chdir('./crosslink_test')


#
# Input Moleclues
#
#'mol1' : 'OC(COC(=O)C=C)COC(=O)C=C'
#'mol2' : 'OCC(COC(=O)C=C)(COC(=O)C=C)COC(=O)C=C'
#'mol3' : 'OCC(COCC(COC(=O)C=C)(COC(=O)C=C)COC(=O)C=C)(COC(=O)C=C)COC(=O)C=C'

l_smiles = {
    'mol1' : 'OCC(COCC(COC(=O)C=C)(COC(=O)C=C)COC(=O)C=C)(COC(=O)C=C)COC(=O)C=C'
    }

l_mols = {
    'mol1' : 100
}

xbox = 70.0
ybox = 70.0
zbox = 70.0

l_reax_atom_type = {
    'reax1' : {
        'mol1' : [12,13,18,19,24,25,30,31,36,37]
    }
}

l_reax_bond_type = {
    'reax1' : {
        'mol1' : [12,18,24,30,36]
    }
}

#
# SMILES to MOL
#
for molname in l_smiles:
    smiles = l_smiles[molname]
    mol_ = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol_)
    AllChem.EmbedMolecule(mol_h)
    mol_3d = Chem.MolToMolBlock(mol_h)
    #
    # Wite MOl File
    #
    f = open(molname+'.mol','w')
    print(mol_3d,file=f)
    f.close()
    #
    # Visualize
    #
    #view(read(molname+'.mol'))
#sys.exit()
#
# Set topology from Molfile
#
# l_topology = {
#    molname : [
#           natom,
#           nbond,
#           bonds: {
#            btype1:
#            btype2:
#            border:
#          }
#       ]
#    }
#
l_topology={}
for molname in l_smiles:
    f = open(molname+'.mol','r')
    lines = f.readlines()
    f.close()
    line = lines[3].split()
    natom = int(line[0])
    nbond = int(line[1])
    topology={}
    topology['natom'] = natom
    topology['nbond'] = nbond
    topology['bonds'] = []
    for line in lines[4+natom:]:
        if 'M  END' in line:
            break
        tmp={}
        line = line.split()
        tmp['btype1'] = int(line[0])
        tmp['btype2'] = int(line[1])
        tmp['border'] = int(line[2])
        topology['bonds'].append(tmp)
    l_topology[molname] = topology
#
# Molecular Packing
#
os.makedirs('packmolfiles',exist_ok=True)
os.chdir('packmolfiles')
g.packmolpath = os.getcwd()
l_moltypes = []
l_molids = []
for molname in l_mols:
    #
    # Wite XYZ File
    #
    read('../'+molname+'.mol').write(molname+'.xyz')
    for i  in range(l_mols[molname]):
        l_moltypes.append(molname)
        xyz = read(molname+'.xyz')
        natom = len(xyz)
        for j in range(natom):
            l_molids.append(i+1)
#print(l_moltypes)
#print(l_molids)    
mk_packmol_random(l_mols,xbox,ybox,zbox)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
#
# make cell
#
cell = read('system.xyz')
lat = cell.get_cell()
lat[0][0] = xbox
lat[1][1] = ybox
lat[2][2] = zbox
cell.set_cell(lat)
cell = sortmol(cell,sort_atom=False)
mk_nvt_input_uff_rigid_scale(cell,0.0005,200,200,xbox/2.0+4.0,20000,20000,400,12345)
os.system('mpirun -np 16 lmp -in lammps.lmp > log_lmp')
tmp = read('result.data',format='lammps-data',sort_by_id=True)
pos = tmp.get_positions(wrap=True)
cell.set_positions(pos)
lat[0][0] = xbox/2.0+4.0
lat[1][1] = ybox/2.0+4.0
lat[2][2] = zbox/2.0+4.0
xbox = lat[0][0]*3.0/4.0
cell.set_cell(lat)
cell.arrays['mol-id'] = l_molids
view(cell)
os.chdir('..')
#
# make full molfile
#
step = 0
f = open('full'+str(step)+'.mol','w')
f.write('system'+'\n')
f.write('  MatStudio         3D                             0'+'\n')
f.write('\n')
f.write('  0  0  0  0  0  0  0  0  0  0999 V3000'+'\n')
f.write('M  V30 BEGIN CTAB'+'\n')

tot_natom = 0
tot_nbond = 0
for moltype in l_moltypes:
    tot_natom = tot_natom + l_topology[moltype]['natom']
    tot_nbond = tot_nbond + l_topology[moltype]['nbond']

f.write('M  V30 COUNTS '+str(tot_natom)+' '+str(tot_nbond)+' 0 0 0'+'\n')
f.write('M  V30 BEGIN ATOM'+'\n')
ista = 0
for moltype in l_moltypes:
    for i in range(ista,ista+l_topology[moltype]['natom']):
        f.write('M  V30 '+str(i+1)+' '+cell[i].symbol+' '+str(cell[i].position[0])+' '+str(cell[i].position[1])+' '+str(cell[i].position[2])+' '+str(i+1)+'\n')
    ista = i+1
f.write('M  V30 END ATOM'+'\n')
f.write('M  V30 BEGIN BOND'+'\n')
ista = 0
natom = 0
bonds = {}
for moltype in l_moltypes:
    j = 0
    for i in range(ista,ista+l_topology[moltype]['nbond']):
        border = l_topology[moltype]['bonds'][j]['border']
        btype1 = l_topology[moltype]['bonds'][j]['btype1']+natom
        btype2 = l_topology[moltype]['bonds'][j]['btype2']+natom
        f.write('M  V30 '+str(i+1)+' '+str(border)+' '+str(btype1)+' '+str(btype2)+'\n')
        bonds[str(i+1)] = [border,btype1,btype2]
        j = j + 1
    ista = i+1
    natom = natom + l_topology[moltype]['natom']
f.write('M  V30 END BOND'+'\n')
f.write('M  V30 BEGIN OBJ3D'+'\n')
f.write('M  V30 END OBJ3D'+'\n')
f.write('M  V30 END CTAB'+'\n')
f.write('M  END'+'\n')
f.close()
#
# Get reaxtype atom-id and bond-id
#
l_reax_atom_id = {}
l_reax_bond_id = {}
for reax in l_reax_atom_type:
    l_reax_atom_id[reax] = []
    #print(list(l_reax_atom_type[reax])[0])
    ista = 0   
    for moltype in l_moltypes:
        j = 0
        for i in range(ista,ista+l_topology[moltype]['natom']):
            if list(l_reax_atom_type[reax])[0] == moltype:
                if j+1 in l_reax_atom_type[reax][moltype]:
                    l_reax_atom_id[reax].append(i+1)
            j = j + 1
        ista = i+1

for reax in l_reax_bond_type:
    ista = 0
    natom = 0
    l_reax_bond_id[reax] = []
    for moltype in l_moltypes:
        j = 0
        for i in range(ista,ista+l_topology[moltype]['nbond']):
            if list(l_reax_bond_type[reax])[0] == moltype:
                if j+1 in l_reax_bond_type[reax][moltype]:
                    btype1 = l_topology[moltype]['bonds'][j]['btype1']+natom
                    btype2 = l_topology[moltype]['bonds'][j]['btype2']+natom
                    l_reax_bond_id[reax].append([i+1,btype1,btype2])
            j = j + 1
        ista = i+1
        natom = natom + l_topology[moltype]['natom']

#for reax in l_reax_atom_id:
#    print(l_reax_atom_id[reax])
#for reax in l_reax_bond_id:
#    print(l_reax_bond_id[reax])

#
# reax_bond:ラジカル化するボンドをランダムに選択
# reax_atom:反応点の原子IDを登録
#
reax_atoms = []
nid = len(l_reax_bond_id['reax1'])
imol = random.randint(0, nid-1)
reax_bonds = l_reax_bond_id['reax1'][imol]
for atom in reax_bonds[1:]:
    reax_atoms.append(atom)
#
#print(reax_atoms)
#
# トポロジーを変える(結合次数2⇒1)
#
bonds[str(reax_bonds[0])][0] = 1
#
# 反応点を抽出し、最短原子間距離の原子idを抽出
#
while True:
    step = step + 1
    if len(reax_atoms) == 1:
        break
    rid1 = random.choice(reax_atoms)
    distances = []
    for reax in l_reax_atom_id:
        for rid2 in l_reax_atom_id[reax]:
            if rid1 != rid2:
                d12_v = cell.get_distance(rid1-1, rid2-1, mic=True, vector=True)
                #if abs(d12_v[0]) > xbox/2.0:
                #    if d12_v[0] < 0.0:
                #        d12_v[0] = d12_v[0] + xbox
                #    elif d12_v[0] > 0.0:
                #        d12_v[0] = d12_v[0] - xbox
                #if abs(d12_v[1]) > ybox/2.0:
                #    if d12_v[1] < 0.0:
                #        d12_v[1] = d12_v[1] + ybox
                #    elif d12_v[1] > 0.0:
                #        d12_v[1] = d12_v[1] - ybox
                #if abs(d12_v[2]) > zbox/2.0:
                #    if d12_v[2] < 0.0:
                #        d12_v[2] = d12_v[2] + zbox
                #    elif d12_v[2] > 0.0:
                #        d12_v[2] = d12_v[2] - zbox
                d12 = math.sqrt(d12_v[0]**2+d12_v[1]**2+d12_v[2]**2)
                distances.append([d12,rid2])
    sorted_d = sorted(distances)
    #while True:
    #    d = random.choice(sorted_d)
    #    if (d[0] < xbox/2.0) and (cell.arrays['mol-id'][rid1-1] == cell.arrays['mol-id'][d[1]-1]):
    #        continue
    #    rid2 = d[1]
    #    d12 = d[0]
    #    print(d12)
    #    break        
    for d in sorted_d:
        if (d[0] < xbox) and (cell.arrays['mol-id'][rid1-1] == cell.arrays['mol-id'][d[1]-1]):
            #print(d[0],d[1])
            continue
        rid2 = d[1]
        d12 = d[0]
        print(d12)
        break
    #print(rid1,rid2,d12)
    #
    # 反応点のトポロジーを変える(結合次数2⇒1)
    #
    for bond in bonds:
        if bonds[bond][0] != 2:
            continue
        if rid2 == bonds[bond][1] or rid2 == bonds[bond][2]:
            bonds[bond][0] = 1
            #
            # 反応点を追加
            #
            for rid_tmp in [bonds[bond][1],bonds[bond][2]]:
                if rid2 != rid_tmp:
                    reax_atoms.append(rid_tmp)
    #
    # 結合をトポロジーに追加
    #
    bonds[str(len(bonds)+1)] = [1,rid1,rid2]
    #
    # 反応に使った原子IDを反応点から削除
    #
    #print(reax_atoms)
    ii = 0
    delid = -1
    for atom in reax_atoms:
        if rid1 == atom:
            delid = ii
        ii = ii + 1
    if delid != -1:
        reax_atoms.pop(delid)
    #ii = 0
    #for atom in reax_atoms:
    #    if rid2 == atom:
    #        delid = ii
    #    ii = ii + 1
    #reax_atoms.pop(delid)

    for reax in l_reax_atom_id:
        ii = 0
        delid = -1
        for atom in l_reax_atom_id[reax]:
            if rid1 == atom:
                delid = ii
            ii = ii + 1
        if delid != -1:
            l_reax_atom_id[reax].pop(delid)

    for reax in l_reax_atom_id:
        ii = 0
        delid = -1
        for atom in l_reax_atom_id[reax]:
            if rid2 == atom:
                delid = ii
            ii = ii + 1
        if delid != -1:
            l_reax_atom_id[reax].pop(delid)

    #print(reax_atoms)
    #
    # 分子IDの更新
    #
    molid1 = cell.arrays['mol-id'][rid1-1]
    molid2 = cell.arrays['mol-id'][rid2-1]
    tmp_molids=[]
    for molid in l_molids:
        if molid == molid2:
            tmp_molids.append(molid1)
        else:
            tmp_molids.append(molid)
    l_molids = tmp_molids
    cell.arrays['mol-id'] = l_molids
    #
    # molファイルの更新
    #
    f = open('full'+str(step)+'.mol','w')
    f.write('system'+'\n')
    f.write('  MatStudio         3D                             0'+'\n')
    f.write('\n')
    f.write('  0  0  0  0  0  0  0  0  0  0999 V3000'+'\n')
    f.write('M  V30 BEGIN CTAB'+'\n')

    tot_nbond = len(bonds)

    f.write('M  V30 COUNTS '+str(tot_natom)+' '+str(tot_nbond)+' 0 0 0'+'\n')
    f.write('M  V30 BEGIN ATOM'+'\n')
    ista = 0
    for moltype in l_moltypes:
        for i in range(ista,ista+l_topology[moltype]['natom']):
            f.write('M  V30 '+str(i+1)+' '+cell[i].symbol+' '+str(cell[i].position[0])+' '+str(cell[i].position[1])+' '+str(cell[i].position[2])+' '+str(i+1)+'\n')
        ista = i+1
    f.write('M  V30 END ATOM'+'\n')
    f.write('M  V30 BEGIN BOND'+'\n')
    i = 0
    for bond in bonds:
        border = bonds[bond][0]
        btype1 = bonds[bond][1]
        btype2 = bonds[bond][2]
        f.write('M  V30 '+str(i+1)+' '+str(border)+' '+str(btype1)+' '+str(btype2)+'\n')
        i = i + 1
    f.write('M  V30 END BOND'+'\n')
    f.write('M  V30 BEGIN OBJ3D'+'\n')
    f.write('M  V30 END OBJ3D'+'\n')
    f.write('M  V30 END CTAB'+'\n')
    f.write('M  END'+'\n')
    f.close()





"""
# 計算セルの長さ(Å)
xbox=100.0
ybox=100.0
zbox=100.0

sigma = 13.55
# モノマーの名前と計算セル中のモノマー分子の数
mols={
    'mol1':16
    }

d_mols = molinit(mols)

ibonds = 1
d_bonds = {}

# モノマーの名前とSMILES式のリスト
smiles={
    'mol1':"OC(COC(=O)C=C)COC(=O)C=C"
    }

reactions={
    'r1':[
        ['$atom:C6','@atom:82'],
        ['$atom:C7','@atom:81L'],
        ['$atom:H18','@atom:85'],
        ['$atom:H19','@atom:85LCH2'],
        ['$atom:H20','@atom:85LCH2']
    ],
    'r2':[
        ['$atom:C12','@atom:82'],
        ['$atom:C13','@atom:81L'],
        ['$atom:H23','@atom:85'],
        ['$atom:H24','@atom:85LCH2'],
        ['$atom:H25','@atom:85LCH2']
    ]
}

reactionids = {
    'mol1':[
        [7,'r1','C6'],
        [8,'r1','C7'],
        [13,'r2','C12'],
        [14,'r2','C13']
    ]
}

# モノマー分子のsmiles式
os.makedirs('rdltfiles',exist_ok=True)
os.chdir('rdltfiles')
g.rdltpath = os.getcwd()
#mollist={}
#bondlist={}
for mol in smiles:
    smi = smiles[mol]
    molname = mol
    # smiles式を読み込む(モノマー)
    smiles2xyz(smi,molname,True,smarts=False,userandom=False)
    monomer = read(molname+'.mol')
    view(monomer)
    cp_rdlt()
    mklt(smi,molname,random=False)
    #mollist0,bondlist0 = rd_lt(rdltpath,molname)
    #for mol0 in mollist0:
    #    mollist[mol0] = mollist0[mol0].copy()
    #for mol0 in bondlist0:
    #    bondlist[mol0] = bondlist0[mol0].copy()
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
g.packmolpath = os.getcwd()
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
g.moltemplatepath = os.getcwd()
shutil.copy(g.packmolpath+'/system.xyz','./')
for mol in smiles:
    molname=mol
    shutil.copy(g.rdltpath+'/'+str(molname)+'.lt','./')
mk_system_lt(d_mols,d_bonds,xbox,ybox,zbox)
os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')
symbols = get_chemical_symbols('system.data')
d_reacatoms = set_reac_atom('./system.data',symbols,reactionids,d_mols)
os.chdir('../')

dis_list = distance_list(g.moltemplatepath+'/system.data',symbols,d_reacatoms,sigma)

a = 59
b = 293

d_mols_tmp = d_mols.copy()
d_mols_tmp = reac_test(a,d_reacatoms,d_mols_tmp,reactions)
d_mols_tmp = reac_test(b,d_reacatoms,d_mols_tmp,reactions)

d_bonds_tmp = d_bonds.copy()
#print(d_reacatoms[a])
#print(d_mols_tmp[d_reacatoms[a]['mol-id']])
#print('$atom:mol'+str(d_reacatoms[a]['mol-id'])+'/'+d_reacatoms[a]['reacatom'])
ibonds_tmp = ibonds
d_bonds_tmp[ibonds_tmp] = [
    '$atom:mol'+str(d_reacatoms[a]['mol-id'])+'[0]/'+d_reacatoms[a]['reacatom'],
    '$atom:mol'+str(d_reacatoms[b]['mol-id'])+'[0]/'+d_reacatoms[b]['reacatom']
]

os.makedirs('tmp',exist_ok=True)
os.chdir('tmp')
shutil.copy(g.packmolpath+'/system.xyz','./')
for imol in d_mols_tmp:
    molname=d_mols_tmp[imol]
    shutil.copy(g.rdltpath+'/'+str(molname)+'.lt','./')
mk_system_lt(d_mols_tmp,d_bonds_tmp,xbox,ybox,zbox)
os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')
symbols = get_chemical_symbols('system.data')
mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,2000,2000,10,10000,400,100,10000,400,200000,300,1.0,'iso',12345,False,qeq=True)
#os.system('mpirun -np 4 lmp -in input_npt.lmp 1> log_lmp 2> err_lmp')
os.chdir('../')
"""

"""
# NPTで圧縮
os.makedirs('npt_compress',exist_ok=True)
os.chdir('npt_compress')
shutil.copy(g.moltemplatepath+'/system.data','./')
shutil.copy(g.moltemplatepath+'/system.in.charges','./')
shutil.copy(g.moltemplatepath+'/system.in.init','./')
shutil.copy(g.moltemplatepath+'/system.in.settings','./')
#### TO DO：NPTの圧縮条件を初期条件へ ####
mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,2000,2000,10,10000,400,100,10000,400,200000,300,1.0,'iso',12345,False,qeq=True)
os.system('mpirun -np 4 lmp -in input_npt.lmp 1> log_lmp 2> err_lmp')
os.chdir('../')
"""



