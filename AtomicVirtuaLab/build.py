def slabgen(cell0,h,k,l,nx,ny,nz,z_lower,z_upper,fname='slab',inv=True):
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.build import make_supercell
    from ase.io import write
    from ase.visualize import view
    cell = AseAtomsAdaptor.get_structure(cell0)
    slabs = SlabGenerator(cell,(h,k,l),nz,10,in_unit_planes=True).get_slabs(repair=True)
    i = 1
    slabs0=[]
    for slab in slabs:
        norm_slab = slab.get_orthogonal_c_slab()
        atoms = AseAtomsAdaptor.get_atoms(norm_slab)
        atoms.wrap()
        atoms_inv = atoms.copy()
        atoms_inv.rotate(180,'x')
        atoms_inv.wrap()
        z=[]
        lat=atoms.get_cell()
        for atom in atoms:
            z.append(atom.position[2])
        zmin=min(z)
        zmax=max(z)
        atoms.translate([0.0,0.0,z_lower-zmin])
        zmax = zmax-zmin
        lat[2][2] = zmax+z_lower+z_upper
        atoms.set_cell(lat)
        atoms = make_supercell(atoms,([nx,0,0],[0,ny,0],[0,0,1]))
        atoms.write(str(fname)+str(i)+'.cif')
        slabs0.append(atoms)
        i = i + 1
        z=[]
        lat=atoms_inv.get_cell()
        for atom in atoms_inv:
            z.append(atom.position[2])
        zmin=min(z)
        zmax=max(z)
        atoms_inv.translate([0.0,0.0,z_lower-zmin])
        zmax = zmax-zmin
        lat[2][2] = zmax+z_lower+z_upper
        atoms_inv.set_cell(lat)
        atoms_inv = make_supercell(atoms_inv,([nx,0,0],[0,ny,0],[0,0,1]))
        if inv:
            atoms_inv.write(str(fname)+str(i)+'.cif')
            slabs0.append(atoms_inv)
            i = i + 1
    return slabs0

def interface_from_slab(substrate_cell0,id0,h_sub,k_sub,l_sub,nx_sub,ny_sub,nz_sub,film_cell0,id1,h_film,k_film,l_film,nx_film,ny_film,nz_film,in_plane_offset=(0.0,0.0),gap=1.6,vacuum_over_film=0.0):
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.interface import Interface
    from ase.visualize import view
    
    #view(substrate_cell0)
    #view(film_cell0)
    
    substrate_cell = AseAtomsAdaptor.get_structure(substrate_cell0)
    substrate_slabs0 = SlabGenerator(substrate_cell,(h_sub,k_sub,l_sub),nz_sub,10,in_unit_planes=True,primitive=False).get_slabs(repair=True)
    film_cell = AseAtomsAdaptor.get_structure(film_cell0)
    film_slabs0 = SlabGenerator(film_cell,(h_film,k_film,l_film),nz_film,10,in_unit_planes=True,primitive=False).get_slabs(repair=True)
    
    #print(len(substrate_slabs0),id0)
    #print(len(film_slabs0),id1)
       
    substrate_slab = substrate_slabs0[id0].get_orthogonal_c_slab()
    film_slab = film_slabs0[id1].get_orthogonal_c_slab()
    
    #substrate_slab = substrate_slabs0[id0].copy()
    #film_slab = film_slabs0[id1].copy()
    
    substrate_slab.make_supercell([nx_sub, ny_sub, 1])
    film_slab.make_supercell([nx_film, ny_film, 1])
    
    #view(AseAtomsAdaptor.get_atoms(substrate_slab))
    #view(AseAtomsAdaptor.get_atoms(film_slab))
          
    interface_slab0 = Interface.from_slabs(substrate_slab, film_slab,in_plane_offset=in_plane_offset,gap=gap,vacuum_over_film=vacuum_over_film)
    #norm_slab = interface_slab0.get_orthogonal_c_slab()
    interface_slab = AseAtomsAdaptor.get_atoms(interface_slab0)
    
    return interface_slab


def sortmol(cell,sort_atom=True):
    from ase import neighborlist
    from ase.build import sort
    from ase.io import write
    from scipy import sparse
    import numpy as np
    if sort_atom:
        cell = sort(cell)
    cutoff = neighborlist.natural_cutoffs(cell)
    neighborList = neighborlist.NeighborList(cutoff, self_interaction=False)
    neighborList.update(cell)
    matrix = neighborList.get_connectivity_matrix()
    n_components, component_list = sparse.csgraph.connected_components(matrix)
    tmp = []
    for l in component_list:
        tmp.append(l+1)
    component_list = np.array(tmp)
    cell.arrays['mol-id'] = component_list
    cell = sort(cell,tags=component_list)
    return cell
    
def sortmol_bond_break(cell,blist):
    from ase import neighborlist
    from ase.build import sort
    from ase.io import write
    from scipy import sparse
    import numpy as np
    cutoff = neighborlist.natural_cutoffs(cell)
    neighborList = neighborlist.NeighborList(cutoff, self_interaction=False, bothways=True)
    neighborList.update(cell)
    matrix = neighborList.get_connectivity_matrix()
    for list in blist:
        matrix[list[0],list[1]] = 0
    n_components, component_list = sparse.csgraph.connected_components(matrix)
    tmp = []
    for l in component_list:
        tmp.append(l+1)
    component_list = np.array(tmp)
    cell.arrays['mol-id'] = component_list
    return component_list,cell