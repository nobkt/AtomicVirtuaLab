from AtomicVirtuaLab.espresso import mk_qe_input_vcrelax, mk_qe_input_scf, mk_qe_input_dos, plot_qe_dos
from AtomicVirtuaLab.io import rd_cif
from ase.io import read
from ase.build import make_supercell
from ase.visualize import view
import AtomicVirtuaLab.globalv as g
import os


g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'
g.fpwo = '/home/A23321P/work/myQE/XSpectra'
g.dos_rusult_dir = '/home/A23321P/work/myQE/XSpectra'

# Si
os.makedirs('./Si',exist_ok=True)
os.chdir('./Si')
#cell = rd_cif(g.cifs+'/Si.cif')
#mk_qe_input_vcrelax(cell,'pbe','paw',level='high',kpts=(4,4,4),ecut='auto')
#cell = read(g.fpwo+'/Si/xanes/qe_scf.pwo')
#mk_qe_input_dos(cell,'pbe','paw',level='high',kpts=(8,8,8),ecut='auto')
#cell = make_supercell(cell,([2,0,0],[0,2,0],[0,0,2]),wrap=True)
#mk_qe_input_scf(cell,'pbe',"paw",level='high',kpts=(2,2,2),ecut='auto')
plot_qe_dos(g.dos_rusult_dir+'/Si/dos',6.2803,nspin=False)
os.chdir('../')


# SiO2
os.makedirs('./SiO2',exist_ok=True)
os.chdir('./SiO2')
#cell = rd_cif(g.cifs+'/SiO2.cif')
#mk_qe_input_vcrelax(cell,'pbe','paw',level='high',kpts=(4,4,4),ecut='auto')
#cell = read(g.fpwo+'/SiO2/xanes/qe_scf.pwo')
#mk_qe_input_dos(cell,'pbe','paw',level='high',kpts=(8,8,8),ecut='auto')
#cell = make_supercell(cell,([2,0,0],[0,2,0],[0,0,2]),wrap=True)
#mk_qe_input_scf(cell,'pbe',"paw",level='high',kpts=(2,2,2),ecut='auto')
plot_qe_dos(g.dos_rusult_dir+'/SiO2/dos',3.6643,nspin=False)
os.chdir('../')
