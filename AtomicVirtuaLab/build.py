def slabgen(cell0,h,k,l,nx,ny,nz,z_lower,z_upper):
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
        atoms.write('slab_'+str(i)+'.cif')
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
        atoms_inv.write('slab_'+str(i)+'.cif')
        slabs0.append(atoms_inv)
        i = i + 1
    return slabs0
