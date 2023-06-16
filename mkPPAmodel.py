mol0 = 'SCi=CC=C(R)C=Ci'

nmer=50

for i in range(nmer):
    j = i + 1
    if j == 1:
        mol = mol0.replace('i',str(j))
    elif j >= 10 and j != nmer:
        mol = mol.replace('R',str(mol0))
        mol = mol.replace('i','%'+str(j))
    elif j < 10 and j != nmer:
        mol = mol.replace('R',str(mol0))
        mol = mol.replace('i',str(j))
    elif j == nmer and j >= 10:
        mol = mol.replace('R',str(mol0))
        mol = mol.replace('i','%'+str(j))
        mol = mol.replace('(R)',"")
    elif j == nmer and j < 10:
        mol = mol.replace('R',str(mol0))
        mol = mol.replace('i',str(j))
        mol = mol.replace('(R)',"")

print(mol)