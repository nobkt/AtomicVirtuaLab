from ase.build import molecule
from ase.visualize import view
from ase.io import read, write
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.build import sortmol
from GenTopo.Coord import PDBobj
from GenTopo.Graph import MolGraph
from GenTopo.GMXTopo import Topo
import os

mollist={
    'H2O':256
}

x_box=100.0
y_box=100.0
z_box=100.0

os.makedirs('test_gromacs',exist_ok=True)
os.chdir('test_gromacs')

mol_ = molecule('H2O')
view(mol_)

os.makedirs('packmol',exist_ok=True)
os.chdir('packmol')

mol_.write('H2O.xyz',format='xyz')
mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cell = read('system.xyz')
cell.set_cell([x_box,y_box,z_box])
view(cell)

os.chdir('..')

cell = sortmol(cell)
cell.new_array('residuenumbers', cell.arrays['mol-id'], int)
cell.set_array('residuenumbers', cell.arrays['mol-id'], int)
cell.write('h2o_box.pdb',format='proteindatabank',write_arrays=True)

mol = PDBobj("h2o_box.pdb")
graph = MolGraph(mol, guessImpropers=True)
gmx = Topo(mol, graph)
gmx.write()
