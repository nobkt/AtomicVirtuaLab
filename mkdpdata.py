from AtomicVirtuaLab.deepmd import qe2dp, get_deepmd_list, wt_deepmd_json
import AtomicVirtuaLab.globalv as g
import sys

g.traindata = "/home/A23321P/work/myPython/AtomicVirtuaLab/deepmd"

#qe2dp(g.traindata,'pwo')
dp_list=get_deepmd_list(g.traindata)
print(dp_list)
#sys.exit()
wt_deepmd_json(g.traindata,dp_list,8.0,1000000,prec='high')
