from ase.data.pubchem import pubchem_atoms_search
from ase.calculators.gaussian import Gaussian
from ase.visualize import view
from ase.io import read
from tangelo import SecondQuantizedMolecule as SQMol
from tangelo.toolboxes.molecular_computation.integral_solver_pyscf import mol_to_pyscf
from tangelo.algorithms.variational import VQESolver, BuiltInAnsatze, SA_VQESolver
from tangelo.algorithms.classical import FCISolver
from tangelo.toolboxes.ansatz_generator.fermionic_operators import spin2_operator
from tangelo.toolboxes.qubit_mappings.mapping_transform import fermion_to_qubit_mapping as f2q_mapping
from pyscf.tools import cubegen
from pyscf import mcscf, gto, scf, lo, dft
from functools import reduce
import numpy
import os
import sys

os.makedirs('./naphthalene_dimer',exist_ok=True)
os.chdir('./naphthalene_dimer')

"""
monomer = pubchem_atoms_search('naphthalene')
#view(monomer)

Gaussian(label='monomer_scf',
         nprocshared=16,
         chk='monomer_scf.chk',
         xc='B3LYP',
         basis='6-31G(d,p)',
         scf='tight,maxcycle=9999',
         opt='calcfc,maxcycles=9999',
         pop='full'
         ).write_input(monomer)

monomer1 = monomer.copy()
monomer2 = monomer.copy()

monomer2.translate([0,0,5])

dimer = monomer1 + monomer2
view(dimer)

Gaussian(label='dimer_scan1',
         nprocshared=16,
         chk='dimer_scan1.chk',
         xc='B3LYP',
         basis='6-31G(d,p)',
         scf='tight,maxcycle=9999',
         opt='ModRedundant,calcfc,maxcycles=9999',
         EmpiricalDispersion='GD3BJ'
         ).write_input(dimer)
"""


# VQE
os.makedirs('VQE',exist_ok=True)
os.chdir('VQE')

"""
# Read Monomer
monomer_tmp = read('/home/A23321P/work/myPySCF/naphthalene_monomer/monomer_opt.log',format='gaussian-out')
#view(monomer_tmp)
monomer = []
for atom in monomer_tmp:
    monomer.append([atom.symbol,(atom.position[0],atom.position[1],atom.position[2])])

# test UHF Natural Orbitals in pySCF
mol_test = gto.M(atom=monomer)
mol_test.basis = '6-31G**'
mol_test.spin = 0
mf = scf.RHF(mol_test)
#mf = scf.UHF(mol_test)
#mf = dft.RKS(mol_test)
#mf = dft.UKS(mol_test)
#mf.xc = 'b3lyp'
mf.max_cycle=2000
mf.kernel()

# Print MO
#print("  #  Energy Occ")
#for i in range(len(mf.mo_energy)):
#    print(f"{i+1:3d}{mf.mo_energy[i]: 9.4f} {int(mf.mo_occ[i])}")

# Output NO to Cube
#for i in range(24,45):
#        cubegen.orbital(mol_test, 'monomer_no_'+str(i)+'.cube', mf.mo_coeff[:, i-1])
#


#mycas = mcscf.CASCI(mf, 2, 2)
mycas = mcscf.CASSCF(mf, 10,(6,4))
mo = mycas.sort_mo([27,31,32,33,34,35,36,37,38,39])
mycas.natorb = True
mycas.kernel(mo)

mycas.verbose = 4
mycas.analyze()

# Print MO
print("  #  Energy Occ")
for i in range(len(mycas.mo_energy)):
    print(f"{i+1:3d}{mycas.mo_energy[i]: 9.4f} {int(mycas.mo_occ[i])}")

# Output NO to Cube
for i in range(30,40):
        cubegen.orbital(mol_test, 'monomer_no_'+str(i)+'.cube', mycas.mo_coeff[:, i-1])

sys.exit()
"""

# Read Gussian Output
mol_tmp = read('/home/A23321P/work/myPySCF/naphthalene_dimer/irc_coords/dimer_irc_f_37.gjf',format='gaussian-in')
#view(mol_tmp)
dimer = []
for atom in mol_tmp:
    dimer.append([atom.symbol,(atom.position[0],atom.position[1],atom.position[2])])

# test UHF Natural Orbitals in pySCF
mol_test = gto.M(atom=dimer)
mol_test.basis = '6-31G**'
mol_test.spin = 0
mf = scf.RHF(mol_test)
#mf = scf.UHF(mol_test)
#mf = dft.UKS(mol_test)
#mf.xc = 'b3lyp'
mf.max_cycle=2000
mf.kernel()

# Print MO
#print("  #  Energy Occ alpha")
#for i in range(len(mf.mo_energy[0])):
#    print(f"{i+1:3d}{mf.mo_energy[0][i]: 9.4f} {int(mf.mo_occ[0][i])}")

#print("  #  Energy Occ beta")
#for i in range(len(mf.mo_energy[1])):
#    print(f"{i+1:3d}{mf.mo_energy[1][i]: 9.4f} {int(mf.mo_occ[1][i])}")

mycas = mcscf.CASCI(mf, 4, 4)
mycas.natorb = True
mycas.kernel()

mycas.verbose = 4
mycas.analyze()

# Print MO
print("  #  Energy Occ")
for i in range(len(mycas.mo_energy)):
    print(f"{i+1:3d}{mycas.mo_energy[i]: 9.4f} {int(mycas.mo_occ[i])}")

# Output NO to Cube
for i in [67,68,69,70]:
        cubegen.orbital(mol_test, 'dimer_no_'+str(i)+'.cube', mycas.mo_coeff[:, i-1])

sys.exit()


# Run RHF
fo = [i for i in range(0,66)]+[j for j in range(70,212)]

mol_dimer = SQMol(dimer, q=0, spin=0, frozen_orbitals=fo, basis="6-31g**")
#mol_dimer_t = SQMol(dimer, q=0, spin=4, frozen_orbitals=fo, basis="3-21g")

# Print MO
print("  #  Energy Occ")
for i in range(len(mol_dimer.mo_energies)):
    print(f"{i+1:3d}{mol_dimer.mo_energies[i]: 9.4f} {int(mol_dimer.mo_occ[i])}")

# Print MO high spin
#print("  #  Energy Occ for high spin")
#for i in range(len(mol_dimer_t.mo_energies)):
#    print(f"{i+1:3d}{mol_dimer_t.mo_energies[i]: 9.4f} {int(mol_dimer_t.mo_occ[i])}")

# Active electrons, Active orbitals
print(f"Number of active electrons: {mol_dimer.n_active_electrons}")
print(f"Number of active orbtials: {mol_dimer.n_active_mos}")

# Output MO to Cube
#for i in [60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]:
#        cubegen.orbital(mol_to_pyscf(mol_dimer,basis = mol_dimer.basis), 'dimer_'+str(i)+'.cube', mol_dimer.mo_coeff[:, i-1])

# Output MO to Cube
#for i in [67,68,69,70]:
#        cubegen.orbital(mol_to_pyscf(mol_dimer_t,basis = mol_dimer_t.basis), 'dimer_t_'+str(i)+'.cube', mol_dimer_t.mo_coeff[:, i-1])

#sys.exit()


# Calculate CASSI
mc = mcscf.CASCI(mol_dimer.mean_field, 4, (2,2))
mc.fcisolver.nroots = 10
mc.kernel()
#mc.verbose = 4
#mc.analyze()


# Dictionary of resources for each algorithm
algorithm_resources = dict()

# Ground state energy calculation with VQE, reference values with FCI
vqe_options = {"molecule": mol_dimer, "ansatz": BuiltInAnsatze.UCCSD}
vqe_solver = VQESolver(vqe_options)
vqe_solver.build()
vqe_energy = vqe_solver.simulate()
print("\n Ground Singlet state")
print(f"VQE energy = {vqe_energy}")
print(f"CASCI energy = {FCISolver(mol_dimer).simulate()}")
algorithm_resources["vqe_ground_state"] = vqe_solver.get_resources()

# Add initial VQE optimal circuit to the deflation circuits list
deflation_circuits = [vqe_solver.optimal_circuit.copy()]

# Calculate first and second excited states by adding optimal circuits to deflation_circuits
for i in range(3):
    vqe_options = {"molecule": mol_dimer, "ansatz": BuiltInAnsatze.UpCCGSD, 
                   "deflation_circuits": deflation_circuits, "deflation_coeff": 0.4}
    vqe_solver = VQESolver(vqe_options)
    vqe_solver.build()
    vqe_energy = vqe_solver.simulate()
    print(f"Excited state #{i+1} \t VQE energy = {vqe_energy}")
    algorithm_resources[f"vqe_deflation_state_{i+1}"] = vqe_solver.get_resources()
    deflation_circuits.append(vqe_solver.optimal_circuit.copy())

# Calculate first and second excited states by adding optimal circuits to deflation_circuits
for i in range(3):
    vqe_options = {"molecule": mol_dimer, "ansatz": BuiltInAnsatze.UpCCGSD, 
               "deflation_circuits": deflation_circuits,
               "deflation_coeff": 0.4, "ref_state": [1, 1, 1, 1, 0, 0, 0, 0]}
    vqe_solver = VQESolver(vqe_options)
    vqe_solver.build()
    vqe_energy = vqe_solver.simulate()
    print(f"Excited state #{i+1} \t VQE energy = {vqe_energy}")
    algorithm_resources[f"vqe_deflation_state_{i+1}"] = vqe_solver.get_resources()
    deflation_circuits.append(vqe_solver.optimal_circuit.copy())
