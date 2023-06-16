def smiles2xyz(smiles,molname,addH):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from ase.io import read, write
    # SMILES式をMOL形式に変換
    m = Chem.MolFromSmiles(smiles)
    # 水素を付加
    if addH:
        m = Chem.AddHs(m)
    #pm = rdDistGeom.ETKDGv2()
    nconf=0
    while True:
        cids = rdDistGeom.EmbedMultipleConfs(m, numConfs=1,maxAttempts=9999,useRandomCoords=True)
        if len(cids) !=0:
            print("len(cids) = "+str(len(cids))+" done!")
            break
        print("len(cids) = "+str(len(cids))+" not searched")
    prop = AllChem.MMFFGetMoleculeProperties(m)
    energy=[]
    for cid in cids:
        mmff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=cid)
        mmff.Minimize()
        energy.append(mmff.CalcEnergy())
    #AllChem.EmbedMolecule(m,p)
    # molファイルとして出力
    Chem.MolToMolFile(m,molname+'.mol')
    # xyzファイルとして出力
    m = read(molname+'.mol')
    m.write(molname+'.xyz')


    