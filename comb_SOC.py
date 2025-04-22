import numpy as np
import sys
from lib_QMout import *

qmout_phm = sys.argv[1]
qmout_ph = sys.argv[2]
qmout_mrci = sys.argv[3]
salida_dir = sys.argv[4]

PH=read_qmout(qmout_ph)
PHM=read_qmout(qmout_phm)
PHPHM=read_qmout(qmout_mrci)

PH_stts,PH_stts_num=getNstates(PH)
PHM_stts,PHM_stts_num=getNstates(PHM)
PHPHM_stts,PHPHM_stts_num=getNstates(PHPHM)

PH_EMat=getE(PH,PH_stts)
PHM_EMat=getE(PHM,PHM_stts)
PHPHM_EMat=getE(PHPHM,PHPHM_stts)

final_mat=map_acoplamientos(PH_EMat,PHM_EMat,PH_stts_num,PHM_stts_num,PHPHM_stts_num,PHPHM_EMat)

dyson_mat=read_Dyson(PHPHM,PHPHM_stts)
xyz=getxyz(PHPHM,PHPHM_stts)
tmx,tmy,tmz=read_tm(PHPHM,PHPHM_stts)

outfile=salida_dir+"/QM.out"
write_output(PHPHM_stts_num,final_mat,dyson_mat,tmx,tmy,tmz,outfile)

