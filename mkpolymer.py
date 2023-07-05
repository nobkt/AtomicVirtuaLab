from AtomicVirtuaLab.io import smiles2molandxyz
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt
from rdkit import Chem
from rdkit.Geometry import Point3D
from ase.io import read
from ase.visualize import view
import os
import codecs

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


