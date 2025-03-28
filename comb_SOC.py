import numpy as np
from lib_QMout import *

PH=read_qmout("PHQM.out")
PHM=read_qmout("PHMQM.out")
PHPHM=read_qmout("PHPHMQM.out")

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

write_output(PHPHM_stts_num,final_mat,dyson_mat,tmx,tmy,tmz)


    

    