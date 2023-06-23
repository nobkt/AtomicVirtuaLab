def smiles2xyz(smiles,molname,addH):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from ase.io import read, write
    import codecs
    # SMILES式をMOL形式に変換
    m = Chem.MolFromSmiles(smiles)
    # 水素を付加
    if addH:
        m = Chem.AddHs(m)
    # 構造最適化
    AllChem.EmbedMolecule(m)
    AllChem.MMFFOptimizeMolecule(m)
    # molファイルとして出力
    Chem.MolToMolFile(m,molname+'.mol')
    # xyzファイルとして出力
    xyz = Chem.MolToXYZBlock(m)
    print(xyz, file=codecs.open(molname+'.xyz', 'w', 'utf-8'))

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






    