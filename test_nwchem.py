from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.nwchem import mk_nwchem_input_opt, mk_nwchem_input_scf
from AtomicVirtuaLab.build import sortmol
import AtomicVirtuaLab.globalv as g
from ase.io import read, write
from ase.visualize import view
from ase.cluster.cubic import FaceCenteredCubic
import os

g.cifdir = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

# Cu-Pc scf
molname='CuPc_3mer'
method='scf'
os.makedirs('./nwchem_test/'+str(molname)+'/'+str(method),exist_ok=True)
os.chdir('./nwchem_test/'+str(molname)+'/'+str(method))

mol_ = read(g.cifdir+'/'+molname+'.xyz')
mol_ = sortmol(mol_,sort_atom=True)
com = mol_.get_center_of_mass()
mol_.translate([-com[0],-com[1],-com[2]])
#view(mol_)

mk_nwchem_input_scf(mol_,'input','B3LYP',basis='6-311G**',chg=0,mult=1)



# Cu-Pc scf end

"""
# nanocluster scf
molname='Pt-cluster'
method='scf'
os.makedirs('./nwchem/'+str(molname)+'/'+str(method),exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/'+str(method))

surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
layers = [6, 9, 5]
lc = 3.9231
ptlayer = FaceCenteredCubic('Pt', surfaces, layers, latticeconstant=lc)
ptlayer.rotate(6, 'x', rotate_cell=True)
ptlayer.rotate(2, 'y', rotate_cell=True)
ptlayer.write('ptlayer.xyz')
mk_nwchem_input_scf(ptlayer,'input','B3LYP',basis='6-31G**',chg=0,mult=1)
view(ptlayer)
x=[]
y=[]
z=[]
for atom in ptlayer:
    x.append(atom.position[0])
    y.append(atom.position[1])
    z.append(atom.position[1])
rx = (max(x)-min(x))/2.0
ry = (max(y)-min(y))/2.0
rz = (max(z)-min(z))/2.0
print(rx,ry,rz,len(ptlayer))
"""

"""
# optimize
molname='PY129'
smiles = 'OC1=CC=CC=C1\\N=C\\C1=C(O)C=CC2=C1C=CC=C2'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY138'
smiles = 'ClC1=C(Cl)C(Cl)=C(Cl)C2=C1C(=O)N(C2=O)C1=CC=CC2=C1NC(C=C2)=C1C(=O)C2=C(C1=O)C(Cl)=C(Cl)C(Cl)=C2Cl'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY139'
smiles = 'O=C1NC(=O)C(=C2NC(C3=C2C=CC=C3)=C2C(=O)NC(=O)NC2=O)C(=O)N1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY150'
smiles = 'OC1=C(\\N=N\\C2=C(O)NC(=O)NC2=O)C(=O)NC(=O)N1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')

molname='PY185'
smiles = 'CNC(=O)C(\\C#N)=C1/NC(C2=C1C=CC=C2)=C1C(=O)NC(=O)NC1=O'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')
"""

"""
molname='PY129_Cu'
smiles = '[#8]-1[Cu]~2[#8]-[#6]-3=[#6](-[#6]=[#7]~2-[#6]-2=[#6]-[#6]=[#6]-[#6]=[#6]-1-2)-[#6]-1=[#6](-[#6]=[#6]-[#6]=[#6]-1)-[#6]=[#6]-3'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True,smarts=True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=2)
os.chdir('../../../')

molname='PY150_Ni'
smiles = 'O=[#6]-1-[#7]-[#6]-2-[#8][Ni]~3[#8]-[#6]-4-[#7]-[#6](=O)-[#7]-[#6](=O)-[#6]-4-[#7]=[#7]~3-[#6]-2-[#6](=O)-[#7]-1'
os.makedirs('./nwchem/'+str(molname)+'/optimize',exist_ok=True)
os.chdir('./nwchem/'+str(molname)+'/optimize')
smiles2xyz(smiles,molname,True,smarts=True)
mol = read(molname+'.xyz')
view(mol)
mk_nwchem_input_opt(mol,molname,'B3LYP',basis='6-31G**',chg=0,mult=1)
os.chdir('../../../')
"""