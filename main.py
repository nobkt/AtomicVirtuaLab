from AtomicVirtuaLab.io import smiles2xyz
from AtomicVirtuaLab.packmol import mk_packmol_random
from AtomicVirtuaLab.moltemplate import cp_rdlt, mklt, mk_system_lt, get_chemical_symbols
from AtomicVirtuaLab.lammps import mk_nvt_input_fr_moltemplate, mk_npt_input_fr_moltemplate,mk_nvtdeform_input_fr_moltemplate
import os

smiles={
    'PR176':'SC1=CC=C(SC2=CC=C(SC3=CC=C(SC4=CC=C(SC5=CC=C(SC6=CC=C(SC7=CC=C(SC8=CC=C(SC9=CC=C(SC%10=CC=C(SC%11=CC=C(SC%12=CC=C(SC%13=CC=C(SC%14=CC=C(SC%15=CC=C(SC%16=CC=C(SC%17=CC=C(SC%18=CC=C(SC%19=CC=C(SC%20=CC=C(SC%21=CC=C(SC%22=CC=C(SC%23=CC=C(SC%24=CC=C(SC%25=CC=C(SC%26=CC=C(SC%27=CC=C(SC%28=CC=C(SC%29=CC=C(SC%30=CC=C(SC%31=CC=C(SC%32=CC=C(SC%33=CC=C(SC%34=CC=C(SC%35=CC=C(SC%36=CC=C(SC%37=CC=C(SC%38=CC=C(SC%39=CC=C(SC%40=CC=C(SC%41=CC=C(SC%42=CC=C(SC%43=CC=C(SC%44=CC=C(SC%45=CC=C(SC%46=CC=C(SC%47=CC=C(SC%48=CC=C(SC%49=CC=C(SC%50=CC=CC=C%50)C=C%49)C=C%48)C=C%47)C=C%46)C=C%45)C=C%44)C=C%43)C=C%42)C=C%41)C=C%40)C=C%39)C=C%38)C=C%37)C=C%36)C=C%35)C=C%34)C=C%33)C=C%32)C=C%31)C=C%30)C=C%29)C=C%28)C=C%27)C=C%26)C=C%25)C=C%24)C=C%23)C=C%22)C=C%21)C=C%20)C=C%19)C=C%18)C=C%17)C=C%16)C=C%15)C=C%14)C=C%13)C=C%12)C=C%11)C=C%10)C=C9)C=C8)C=C7)C=C6)C=C5)C=C4)C=C3)C=C2)C=C1'
}

mollist={
    'PR176':200
}

x_box=1000.0
y_box=1000.0
z_box=1000.0



os.makedirs('./packmol',exist_ok=True)
os.chdir('./packmol')

for smi in smiles:
    smiles2xyz(smiles[smi],smi,True)

mk_packmol_random(mollist,x_box,y_box,z_box)
os.system('packmol < packmol.inp 1> log_packmol 2> err_packmol')

cp_rdlt()

for smi in smiles:
    mklt(smiles[smi],smi)

mk_system_lt(mollist,x_box,y_box,z_box)

os.system('moltemplate.sh -xyz system.xyz -nocheck system.lt')
os.system('cleanup_moltemplate.sh')

symbols = get_chemical_symbols('system.data')

#mk_nvt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,12345,False)
mk_npt_input_fr_moltemplate(symbols,True,1.0,200,200,20000,300,0,'iso',12345,False)
mk_nvtdeform_input_fr_moltemplate(symbols,False,1.0,200,200,20000,300,1.0e10*1.0e-15,200,12345,True)

os.chdir('../')

