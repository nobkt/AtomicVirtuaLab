from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mklt_fr_mol, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.lammps import mk_npt_input_fr_moltemplate, mk_npt_compress_input_fr_moltemplate
import AtomicVirtuaLab.globalv as g
from ase.visualize import view
from ase.io import read
from ase.geometry.analysis import Analysis
import pandas as pd
import sys
from rdkit import Chem
from rdkit.Chem import Draw, rdDistGeom
import os
import shutil
import sys

g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'



# モノマー溶液MD
nmol = 128

os.makedirs('./uv_risin_md',exist_ok=True)
os.chdir('./uv_risin_md')

#mols = Chem.SDMolSupplier(g.cifs+'/uv_risin_molecules.sdf')
#mols = Chem.SDMolSupplier(g.cifs+'/uv_risin_molecules2.sdf')
#mols = Chem.SDMolSupplier(g.cifs+'/uv_risin_new20231124.sdf')
mols = Chem.SDMolSupplier(g.cifs+'/uv_risin_new20231124_2.sdf')

i = 1
for mol in mols:
    if i==7:
        rnd = True
    else:
        rnd = True
    os.makedirs('mol'+str(i),exist_ok=True)
    os.chdir('mol'+str(i))
    smiles=Chem.MolToSmiles(mol)
    print('mol'+str(i)+' : ',smiles)
    cp_rdlt()
    mklt(smiles,'mol'+str(i),random=rnd)
    #smiles2xyz(smiles,'mol'+str(i),True,smarts=False,userandom=False)
    Draw.MolToFile(mol,'mol'+str(i)+'.png')
    mollist={
        'mol'+str(i):nmol
    }
    x_box=300.0
    y_box=300.0
    z_box=300.0
    os.makedirs('./packmol',exist_ok=True)
    os.chdir('./packmol')
    smiles2xyz(smiles,'mol'+str(i),True,smarts=False,userandom=rnd)
    view(read('mol'+str(i)+'.xyz'))
    mk_packmol_random(mollist,x_box,y_box,z_box)
    os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')
    #tmp = read('system.xyz')
    #tmp.set_cell([x_box,y_box,z_box])
    #view(tmp)
    #tmp.write('test.cif')
    packmolpath = os.getcwd()
    os.chdir('../')
    os.makedirs('moltemplate',exist_ok=True)
    os.chdir('moltemplate')
    shutil.copy('../packmol/system.xyz','./')
    shutil.copy('../mol'+str(i)+'.lt','./')
    mk_system_lt(mollist,x_box,y_box,z_box)
    os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
    os.system('cleanup_moltemplate.sh')
    symbols = get_chemical_symbols('system.data')
    os.chdir('../')
    os.makedirs('npt_p100',exist_ok=True)
    os.chdir('npt_p100')
    shutil.copy('../moltemplate/system.data','./')
    shutil.copy('../moltemplate/system.in.charges','./')
    shutil.copy('../moltemplate/system.in.init','./')
    shutil.copy('../moltemplate/system.in.settings','./')
    #mk_npt_input_fr_moltemplate(symbols,True,0.5,200,200,200000,400,100,'iso',12345,False,qeq=True)
    f = open('system.data','r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if 'atom types' in line:
            line = line.split()
            ntype = int(line[0])
    rdfpairs=[]
    for ipair in range(ntype-1):
        for jpair in range(ipair,ntype):
            rdfpairs.append([ipair+1,jpair+1])
    mk_npt_compress_input_fr_moltemplate(symbols,True,0.5,2000,2000,10,10000,400,100,10000,400,2000000,300,1.0,'iso',12345,False,qeq=True,rdfpairs=rdfpairs)
    os.chdir('../')
    i+=1
    os.chdir('../')
# モノマー溶液MD 終了


"""
# RDF解析
cell = read('./traj_npt.lammpstrj',format='lammps-dump-text',index=':')
intra=True
inter=True
neq=20
type1=4
type2=4

cell0=[]
for atom in cell:
  atom.wrap()
  cell0.append(atom)

i = 0
nn=len(cell0)
dn = int(nn/neq)
for n in range(0,nn,dn):
  atoms = cell0[n]
  i = i + 1
  if intra:
      distribution1, distance = Analysis(atoms).get_rdf_intramol(rmax=15., nbins=100, return_dists=True, types=[type1,type2])[0]
  if inter:
      distribution2, distance = Analysis(atoms).get_rdf_intermol(rmax=15., nbins=100, return_dists=True, types=[type1,type2])[0]
  if i == 1:
    if intra:
        intra_distribution = distribution1
    if inter:
        inter_distribution = distribution2
  else:
    if intra:
        intra_distribution = [x+y for (x,y) in zip(intra_distribution,distribution1)]
    if inter:
        inter_distribution = [x+y for (x,y) in zip(inter_distribution,distribution2)]

if intra:
    intra_distribution = [x/float(i) for x in intra_distribution]
if inter:
    inter_distribution = [x/float(i) for x in inter_distribution]

if intra:
    intra_df = pd.DataFrame({'distance':distance,'rdf':intra_distribution})
    intra_df.to_csv('intra_rdf_'+str(type1)+'_'+str(type2)+'.csv')
if inter:
    inter_df = pd.DataFrame({'distance':distance,'rdf':inter_distribution})
    inter_df.to_csv('inter_rdf_'+str(type1)+'_'+str(type2)+'.csv')
# RDF解析 終了
"""


