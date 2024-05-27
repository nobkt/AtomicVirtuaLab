def smiles2xyz(smiles,molname,addH,smarts=False,userandom=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from ase.io import read, write
    import codecs
    # SMILES式をMOL形式に変換
    if smarts:
        m = Chem.MolFromSmarts(smiles)
    else:
        m = Chem.MolFromSmiles(smiles)
    # 水素を付加
    if addH:
        m = Chem.AddHs(m)
    # 構造最適化
    if userandom:
        AllChem.EmbedMolecule(m,useRandomCoords=True,randomSeed = 12345)
        AllChem.MMFFOptimizeMolecule(m)
    else:
        AllChem.EmbedMolecule(m)
        AllChem.MMFFOptimizeMolecule(m)
    # molファイルとして出力
    Chem.MolToMolFile(m,molname+'.mol')
    # xyzファイルとして出力
    xyz = Chem.MolToXYZBlock(m)
    print(xyz, file=codecs.open(molname+'.xyz', 'w', 'utf-8'))

def smiles2molandxyz(smiles,molname,addH,smarts=False,userandom=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from ase.io import read, write
    import codecs
    # SMILES式をMOL形式に変換
    if smarts:
        m = Chem.MolFromSmarts(smiles)
    else:
        m = Chem.MolFromSmiles(smiles)
    # 水素を付加
    if addH:
        m = Chem.AddHs(m)
    # 構造最適化
    if userandom:
        AllChem.EmbedMolecule(m,useRandomCoords=True,randomSeed = 12345)
        AllChem.MMFFOptimizeMolecule(m)
    else:
        AllChem.EmbedMolecule(m)
        AllChem.MMFFOptimizeMolecule(m)
    # molファイルとして出力
    Chem.MolToMolFile(m,molname+'.mol')
    # xyzファイルとして出力
    xyz = Chem.MolToXYZBlock(m)
    print(xyz, file=codecs.open(molname+'.xyz', 'w', 'utf-8'))
    return m

def smiles2molandpdb(smiles,molname,addH,smarts=False):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from ase.io import read, write
    import codecs
    # SMILES式をMOL形式に変換
    if smarts:
        m = Chem.MolFromSmarts(smiles)
    else:
        m = Chem.MolFromSmiles(smiles)
    # 水素を付加
    if addH:
        m = Chem.AddHs(m)
    # 構造最適化
    AllChem.EmbedMolecule(m)
    AllChem.MMFFOptimizeMolecule(m)
    # molファイルとして出力
    Chem.MolToMolFile(m,molname+'.mol')
    # pdbファイルとして出力
    writer = Chem.PDBWriter(molname+'.pdb')
    writer.write(m)
    return m

def rd_cif(fcif,primitive_cell=False):
    from ase.io import read
    if primitive_cell:
        cell = read(fcif,primitive_cell=True, subtrans_included=False)
    else:
        cell = read(fcif)
    return cell

def rd_contcar(fcontcar):
    from ase.io import read
    cell = read(fcontcar)
    return cell

def cell2atomlist(cell):
    import re
    from ase.build import sort
    cell0 = sort(cell)
    symbols0 = cell0.get_chemical_symbols()
    symbols=[]
    for symbol in symbols0:
        if symbol in symbols:
            continue
        else:
            symbols.append(symbol)
    #symbols = re.split('[0-9]',formula)
    #symbols = [l for l in symbols if l != '']
    return symbols

def mk_lammpsdata(cell,charge,force_skew=False,mol=False):
    from ase.io import write
    if charge:
        cell.write('lammps.data',format='lammps-data',force_skew=force_skew,atom_style='charge')
    elif mol:
        cell.write('lammps.data',format='lammps-data',force_skew=force_skew,atom_style='full')
    else:
        cell.write('lammps.data',format='lammps-data',force_skew=force_skew)

def rd_lammpsdata(cell,fdata,charge):
    from ase.io import read
    symbols = cell.get_chemical_symbols()
    if charge:
        cell = read(fdata,format='lammps-data',style='charge',sort_by_id=True)
    else:
        cell = read(fdata,format='lammps-data',style='atomic',sort_by_id=True)
    cell.set_chemical_symbols(symbols)
    return cell

def rd_lammpsdata_init(Z_of_type,fdata,charge):
    from ase.io import read
    if charge:
        cell = read(fdata,format='lammps-data',style='full',sort_by_id=True,Z_of_type=Z_of_type)
    else:
        cell = read(fdata,format='lammps-data',style='atomic',sort_by_id=True,Z_of_type=Z_of_type)
    return cell

   