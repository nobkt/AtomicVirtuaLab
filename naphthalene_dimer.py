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
from tangelo.toolboxes.operators import FermionOperator
from pyscf.tools import cubegen
from pyscf import mcscf, gto, scf, lo, dft
from openfermion.utils import hermitian_conjugated as hc
from scipy.linalg import eigh
from functools import reduce
import numpy
import numpy as np
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
mol_test.basis = '3-21G'
mol_test.spin = 2
mol_test.build()
mf = scf.RHF(mol_test)
#mf = scf.UHF(mol_test)
#mf = dft.RKS(mol_test)
#mf = dft.UKS(mol_test)
#mf.xc = 'b3lyp'
mf.max_cycle=2000
mf.kernel()

# Print MO
print("  #  Energy Occ")
for i in range(len(mf.mo_energy)):
    print(f"{i+1:3d}{mf.mo_energy[i]: 9.4f} {int(mf.mo_occ[i])}")

# Output MO to Cube
#for i in range(20,50):
#        cubegen.orbital(mol_test, 'monomer_no_'+str(i)+'.cube', mf.mo_coeff[:, i-1])
#sys.exit()

#mycas = mcscf.CASCI(mf, 2, 2)
mycas = mcscf.CASCI(mf, 10,10)
mo = mycas.sort_mo([27,31,32,33,34,35,36,37,40,47])
mycas.natorb = True
mycas.fcisolver.nroots = 20
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


"""
# test UHF Natural Orbitals in pySCF
mol_test = gto.M(atom=dimer)
mol_test.basis = '3-21g'
mol_test.spin = 0
mol_test.build()
mf = scf.RHF(mol_test)
#mf = scf.UHF(mol_test)
#mf = dft.UKS(mol_test)
#mf.xc = 'b3lyp'
mf.max_cycle=2000
mf.kernel()

# Print MO
#print("  #  Energy Occ alpha")
#for i in range(len(mf.mo_energy)):
#    print(f"{i+1:3d}{mf.mo_energy[i]: 9.4f} {int(mf.mo_occ[i])}")

# Output MO to Cube
#for i in [65,66,67,68,69,70,71,72]:
#        cubegen.orbital(mol_test, 'dimer_no_'+str(i)+'.cube', mf.mo_coeff[:, i-1])
#sys.exit()

mycas = mcscf.CASCI(mf, 4, 4)
mycas.fcisolver.nroots = 30
mycas.natorb = True
mycas.kernel()

mycas.verbose = 4
mycas.analyze()

# Print MO
print("  #  Energy Occ")
for i in range(len(mycas.mo_energy)):
    print(f"{i+1:3d}{mycas.mo_energy[i]: 9.4f} {int(mycas.mo_occ[i])}")

# Output NO to Cube
for i in [65,66,67,68,69,70,71,72]:
        cubegen.orbital(mol_test, 'dimer_no_'+str(i)+'.cube', mycas.mo_coeff[:, i-1])

sys.exit()
"""

# Run RHF
fo = [i for i in range(0,66)]+[j for j in range(70,212)]

mol_dimer = SQMol(dimer, q=0, spin=0, frozen_orbitals=fo, basis="3-21g")
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
#        cubegen.orbital(mol_to_pyscf(mol_dimer,basis = mol_dimer.basis), 'dimer_'+str(i)+'.cube', mol_dimer.mo_coeff[:, i-1])

#sys.exit()

# Calculate CASSI
#mc = mcscf.CASCI(mol_dimer.mean_field, 8,8)
#mc.natorb = True
#mc.fcisolver.nroots = 30
#mc.kernel()
#mc.verbose = 4
#mc.analyze()
#sys.exit()


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
"""

# Add initial VQE optimal circuit to the deflation circuits list
deflation_circuits = [vqe_solver.optimal_circuit.copy()]

# Calculate first and second excited states by adding optimal circuits to deflation_circuits
for i in range(8):
    vqe_options = {"molecule": mol_dimer, "ansatz": BuiltInAnsatze.UpCCGSD, 
                   "deflation_circuits": deflation_circuits, "deflation_coeff": 0.4}
    vqe_solver = VQESolver(vqe_options)
    vqe_solver.build()
    vqe_energy = vqe_solver.simulate()
    print(f"Excited state #{i+1} \t VQE energy = {vqe_energy}")
    algorithm_resources[f"vqe_deflation_state_{i+1}"] = vqe_solver.get_resources()
    deflation_circuits.append(vqe_solver.optimal_circuit.copy())

# Calculate first and second excited states by adding optimal circuits to deflation_circuits
vqe_options = {"molecule": mol_dimer, "ansatz": BuiltInAnsatze.UpCCGSD, 
           "deflation_circuits": deflation_circuits,
           "deflation_coeff": 0.4, "ref_state": [1,0,0,0,1,0,0,0]}
vqe_solver = VQESolver(vqe_options)
vqe_solver.build()
vqe_energy = vqe_solver.simulate()
print(f"Excited state #{i+1} \t VQE energy = {vqe_energy}")
algorithm_resources[f"vqe_deflation_state_{9}"] = vqe_solver.get_resources()
deflation_circuits.append(vqe_solver.optimal_circuit.copy())

sys.exit()
"""

"""
#State-Averaged VQE
algorithm_resources = dict()

vqe_options = {"molecule": mol_dimer, "ref_states": [[1,0,0,1,1,0,0,1]],
               "weights": [1.0], "penalty_terms": None,
               "qubit_mapping": "jw", "ansatz": BuiltInAnsatze.UpCCGSD,
               }
vqe_solver = SA_VQESolver(vqe_options)
vqe_solver.build()
enernew = vqe_solver.simulate()
for i, energy in enumerate(vqe_solver.state_energies):
    print(f"Singlet State {i} has energy {energy}")

algorithm_resources["sa_vqe"] = vqe_solver.get_resources()

s2op = f2q_mapping(spin2_operator(4), "jw")
for i in range(1):
    print(f"State {i} has S^2 = {vqe_solver.backend.get_expectation_value(s2op, vqe_solver.reference_circuits[i]+vqe_solver.optimal_circuit)}")
sys.exit()

#Multistate, contracted VQE (MC-VQE)
# Generate individual statevectors
ref_svs = list()
for circuit in vqe_solver.reference_circuits:
    _, sv = vqe_solver.backend.simulate(circuit, return_statevector=True)
    ref_svs.append(sv)

# Generate Equation (2) using equation (4) and (5) of arXiv:1901.01234
h_theta_theta = np.zeros((2,2))
for i, sv1 in enumerate(ref_svs):
    for j, sv2 in enumerate(ref_svs):
        if i != j:
            sv_plus = (sv1 + sv2)/np.sqrt(2)
            sv_minus = (sv1 - sv2)/np.sqrt(2)
            exp_plus = vqe_solver.backend.get_expectation_value(vqe_solver.qubit_hamiltonian, vqe_solver.optimal_circuit, initial_statevector=sv_plus)
            exp_minus = vqe_solver.backend.get_expectation_value(vqe_solver.qubit_hamiltonian, vqe_solver.optimal_circuit, initial_statevector=sv_minus)
            h_theta_theta[i, j] = (exp_plus-exp_minus)/2
        else:
            h_theta_theta[i, j] = vqe_solver.state_energies[i]

e, _ = np.linalg.eigh(h_theta_theta)
for i, energy in enumerate(e):
    print(f"Singlet State {i} \t MC-VQE energy = {energy}")
"""

# Quantum Subspace Expansion
# Generate all single excitations as qubit operators
op_list = list()
#for i in range(2):
#    for j in range(i+1, 2):
#        op_list += [f2q_mapping(FermionOperator(((2*i, 1), (2*j, 0))), "jw")] #spin-up transition
#        op_list += [f2q_mapping(FermionOperator(((2*i+1, 1), (2*j+1, 0))), "jw")] #spin-down transition
#        op_list += [f2q_mapping(FermionOperator(((2*i+1, 1), (2*j, 0))), "jw")] #spin-up to spin-down
#        op_list += [f2q_mapping(FermionOperator(((2*i, 1), (2*j+1, 0))), "jw")] #spin-down to spin-up

# single exitations
for i in range(4):
    for j in range(4):
        op_list += [f2q_mapping(FermionOperator(((i, 1), (j+4, 0))), "jw")]
        print(i,j+4)

# double exitations
for i in range(4):
    for j in range(i+1,4):
        for k in range(4):
            for l in range(k+1,4):
                op_list += [f2q_mapping(FermionOperator(((i, 1), (j, 1), (k+4, 0), (l+4, 0))), "jw")]
                print(i,j,k+4,l+4)


# triple exitations
for i in range(4):
    for j in range(i+1,4):
        for k in range(j+1,4):
            for l in range(4):
                for m in range(l+1,4):
                    for n in range(m+1,4):
                        op_list += [f2q_mapping(FermionOperator(((i, 1), (j, 1), (k, 1), (l+4, 0), (m+4, 0), (n+4, 0))), "jw")]
                        print(i,j,k,l+4,m+4,n+4)

# 4 exitations
for i in range(4):
    for j in range(i+1,4):
        for k in range(j+1,4):
            for l in range(k+1,4):
                for m in range(4):
                    for n in range(m+1,4):
                        for o in range(n+1,4):
                            for p in range(o+1,4):
                                op_list += [f2q_mapping(FermionOperator(((i, 1), (j, 1), (k, 1), (l, 1), (m+4, 0), (n+4, 0), (o+4, 0), (p+4, 0))), "jw")]
                                print(i,j,k,l,m+4,n+4,o+4,p+4)

# Compute F and S matrices.
size_mat = len(op_list)
h = np.zeros((size_mat, size_mat))
s = np.zeros((size_mat, size_mat))
state_circuit = vqe_solver.optimal_circuit
for i, op1 in enumerate(op_list):
    for j, op2 in enumerate(op_list):
        h[i, j] = np.real(vqe_solver.backend.get_expectation_value(hc(op1)*vqe_solver.qubit_hamiltonian*op2, state_circuit))
        s[i, j] = np.real(vqe_solver.backend.get_expectation_value(hc(op1)*op2, state_circuit))
#
label = "quantum_subspace_expansion"
algorithm_resources[label] = vqe_solver.get_resources()
algorithm_resources[label]["n_post_terms"] = len(op_list)**2*algorithm_resources[label]["qubit_hamiltonian_terms"]

# Solve FU = SUE
e, v = eigh(h,s)
print(f"Quantum Subspace Expansion energies: \n {e}")
#print(v)
