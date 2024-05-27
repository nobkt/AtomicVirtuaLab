def get_density(cell):
    from ase.units import mol
    tot_mass = sum(cell.get_masses())
    volume = cell.get_volume()
    density = (tot_mass/mol)/(volume*1.0e-24)
    return density

def scale_cell(cell,scale_factor):
    lat = cell.get_cell_lengths_and_angles()
    lat[0] = lat[0] * scale_factor
    lat[1] = lat[1] * scale_factor
    lat[2] = lat[2] * scale_factor
    cell0 = cell.copy()
    cell0.set_cell(lat,scale_atoms=True)
    return cell0

def add_displacement(cell,dr):
    import random
    positions0 = cell.get_positions()
    positions=[]
    for pos in positions0:
        x = random.gauss(pos[0],dr)
        y = random.gauss(pos[1],dr)
        z = random.gauss(pos[2],dr)
        positions.append([x,y,z])
    cell0 = cell.copy()
    cell0.set_positions(positions)
    cell0.wrap()
    return cell0

def deform_cell(cell,da,dtheta):
    import random
    lat = cell.get_cell_lengths_and_angles()
    a = random.gauss(lat[0],da)
    b = random.gauss(lat[1],da)
    c = random.gauss(lat[2],da)
    alpha = random.gauss(lat[3],dtheta)
    beta = random.gauss(lat[4],dtheta)
    gamma = random.gauss(lat[5],dtheta)
    cell0 = cell.copy()
    cell0.set_cell([a,b,c,alpha,beta,gamma],scale_atoms=True)
    return cell0

