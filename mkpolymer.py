from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt
from AtomicVirtuaLab.build import sortmol_bond_break
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Geometry import Point3D
from ase.io import read, write
from ase.visualize import view
import os
import codecs
import numpy as np
import sys



os.makedirs('./polymer_build',exist_ok=True)
os.chdir('./polymer_build')

#smiles_3mer = 'SC1=CC=C(SC2=CC=C(SC3=CC=CC=C3)C=C2)C=C1'
smiles_3mer = 'C[Si](C)(O)O[Si](O)(OC1=CC=CC=C1)O[Si](O)(O)OCCOCC1CO1'

cp_rdlt()
#mklt(smiles_3mer,'3mer',random=True)

mol_3mer = smiles2xyz(smiles_3mer,'3mer',True,smarts=False,userandom=False)
sys.exit()
cell = read('./3mer.mol')
#view(cell)

blist=[
    [4,5],
    [5,4],
    [9,10],
    [10,9]
]

component_list,cell = sortmol_bond_break(cell,blist)

#print(component_list)
cell.set_tags(component_list)

f = open('./3mer.lt','r')
lines = f.readlines()
f.close()
f = open('./test.lt','w')
for line in lines:
    f.write(line)
f.close()

mol1 = cell.copy()
del mol1[[atom.index for atom in mol1 if atom.tag != 1]]
mol1.write('mol1.xyz',format='xyz')
#view(mol1)

mol2 = cell.copy()
del mol2[[atom.index for atom in mol2 if atom.tag != 2]]
mol2.write('mol2.xyz',format='xyz')
#view(mol2)

mol3 = cell.copy()
del mol3[[atom.index for atom in mol3 if atom.tag != 3]]
mol3.write('mol3.xyz',format='xyz')
#view(mol3)


"""
rxn = AllChem.ReactionFromSmarts("[Rb][*:1].[Ra][*:2]>>[*:1][*:2]")
rxn0 = AllChem.ReactionFromSmarts("[Ra][*:1].[Rb][*:2]>>[*:1][*:2]")
s1 = "c1([Ra])ccc(S[Rb])cc1"
m1 = Chem.MolFromSmiles(s1)
m1 = Chem.AddHs(m1)
AllChem.EmbedMolecule(m1)
s2= "c1([Ra])ccc(S[Rb])cc1"
m2 = Chem.MolFromSmiles(s2)
m2 = Chem.AddHs(m2)
AllChem.EmbedMolecule(m2)
s0 = "c1ccc(S[Rb])cc1"
m0 = Chem.MolFromSmiles(s0)
m0 = Chem.AddHs(m0)
AllChem.EmbedMolecule(m0)
s4 = "c1([Ra])ccc(S[H])cc1"
m4 = Chem.MolFromSmiles(s4)
m4 = Chem.AddHs(m4)
AllChem.EmbedMolecule(m4)
results = rxn.RunReactants( [m1, m2] )
for i in range(50):
    for products in results:
        for mol in products:
            #mol.RemoveAllConformers()
            #Chem.rdCoordGen.AddCoords(mol)
            results = rxn.RunReactants( [mol, m2] )
for products in results:
    for mol in products:
        #mol.RemoveAllConformers()
        #Chem.rdCoordGen.AddCoords(mol)
        results = rxn0.RunReactants( [mol, m0] )
for products in results:
    for mol in products:
        #mol.RemoveAllConformers()
        #Chem.rdCoordGen.AddCoords(mol)
        results = rxn.RunReactants( [mol, m4] )
for products in results:
    for mol in products:
        AllChem.EmbedMolecule(mol,useRandomCoords=True)
        #AllChem.MMFFOptimizeMolecule(mol)
        Chem.MolToMolFile(mol,'test.mol')
        cp_rdlt()
        mklt(s1,'3mer')

      
mol = read('./test.mol')
view(mol)

#print(cids)
"""


"""
smiles = 'OCC1OC(OC2C(O)C(O)COC2CO)C(O)C(O)C1O'

os.makedirs('./polymer_build',exist_ok=True)
os.chdir('./polymer_build')

cp_rdlt()
mklt(smiles,'PPS1mer')

mol = smiles2molandxyz(smiles,'PPS1mer',True)

mol_tmp = Chem.Mol(mol)

view(read('PPS1mer.mol'))

id_tail = 10
id_head = 7

position_head = mol.GetConformer().GetAtomPosition(id_head)
position_tail = mol.GetConformer().GetAtomPosition(id_tail)

connect_head = mol_tmp.GetAtomWithIdx(id_head).GetNeighbors()
for atom in connect_head:
    id_connect_head = atom.GetIdx()
    position_connect_head = mol_tmp.GetConformer().GetAtomPosition(id_connect_head)

connect_tail = mol_tmp.GetAtomWithIdx(id_tail).GetNeighbors()
for atom in connect_tail:
    id_connect_tail = atom.GetIdx()
    position_connect_tail = mol_tmp.GetConformer().GetAtomPosition(id_connect_tail)

print(id_connect_head,position_connect_head.x,position_connect_head.y,position_connect_head.z)
print(id_connect_tail,position_connect_tail.x,position_connect_tail.y,position_connect_tail.z)

t_vec = [position_head.x - position_connect_tail.x,
         position_head.y - position_connect_tail.y,
         position_head.z - position_connect_tail.z]

for i in range(len(mol_tmp.GetAtoms())):
    p = mol_tmp.GetConformer().GetAtomPosition(i)
    mol_tmp.GetConformer().SetAtomPosition(i, Point3D(p[0]+t_vec[0], p[1]+t_vec[1], p[2]+t_vec[2]))




#Chem.MolToMolFile(test_mol,'test'+'.mol')
#xyz = Chem.MolToXYZBlock(test_mol)
#print(xyz, file=codecs.open('test'+'.xyz', 'w', 'utf-8'))
#print(Chem.MolToMolBlock(mol))
#print(Chem.MolToMolBlock(mol2))
#view(read('test.mol'))
#view(read('test.xyz'))
"""

