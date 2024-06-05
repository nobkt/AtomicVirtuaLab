from ase.io import read
from ase.visualize import view
import os

for U in range(10):
        os.chdir('U='+str(U))
        cell = read('qe_vc-relax.pwo')
        #if U == 0:
        #       view(cell)
        #lat = cell.get_cell_lengths_and_angles()
        #print(lat)
        d1=[]
        d2=[]
        for i in [20,32,16,18,30,34]:
                d = cell.get_distance(4, i, mic=False, vector=False)
                d2.append(d)
                if i == 20 or i== 32:
                        d1.append(d)
        #       print(U,4,i,d)
        print(U,d1)
        lave = sum(d2)/len(d2)
        print(U,lave)
        d0=0.0
        for d in d2:
                d0 = d0 + abs(d-lave)/lave
        print(U,d0/len(d2))
        d1=[]
        d2=[]
        for i in [21,33,17,19,31,35]:
                d = cell.get_distance(5, i, mic=False, vector=False)
                d2.append(d)
                if i == 21 or i == 33:
                        d1.append(d)
        #       print(U,5,i,d)
        print(U,d1)
        lave = sum(d2)/len(d2)
        print(U,lave)
        d0=0.0
        for d in d2:
                d0 = d0 + abs(d-lave)/lave
        print(U,d0/len(d2))
        os.chdir('..')
        
