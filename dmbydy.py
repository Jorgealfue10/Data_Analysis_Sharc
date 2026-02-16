import sys
from lib_QMout import * 

inputfile="./QM.out"
out_file="./QM_test.out"

data=read_qmout(inputfile)
total_stts,stts=getNstates(data)
e_mat=getE(data,total_stts)
tmx,tmy,tmz=read_tm(data,total_stts)
dyson=read_Dyson(data,total_stts)
e_mat=np.array(e_mat)
write_NoSOC_QMout(data,stts,e_mat,tmx,tmy,tmz,dyson,out_file)
