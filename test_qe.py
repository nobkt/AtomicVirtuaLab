from AtomicVirtuaLab.espresso import mk_qe_input_vcrelax, mk_qe_input_dos, mk_qe_input_band
from AtomicVirtuaLab.io import rd_cif
import AtomicVirtuaLab.globalv as g
import os

g.qepot = '/home/A23321P/work/myPython/AtomicVirtuaLab/qe_pseudo'
g.cifs = '/home/A23321P/work/myPython/AtomicVirtuaLab/cifs'

cell = rd_cif(g.cifs+'/'+'CuI_zincblende.cif',primitive_cell=True)

level='high'
os.makedirs('./qe/vc-relax/'+str(level),exist_ok=True)
os.chdir('./qe/vc-relax/'+str(level))
mk_qe_input_vcrelax(cell,'pbe','paw',kpts=(4,4,4),ecut='auto',level=str(level))
os.chdir('../../..')
os.makedirs('./qe/dos/'+str(level),exist_ok=True)
os.chdir('./qe/dos/'+str(level))
mk_qe_input_dos(cell,'pbe','paw',level=str(level),kpts=(10,10,10),ecut='auto')
os.chdir('../../..')
os.makedirs('./qe/band/'+str(level),exist_ok=True)
os.chdir('./qe/band/'+str(level))
mk_qe_input_band(cell,'pbe','paw',level=str(level),ecut='auto')
os.chdir('../../..')