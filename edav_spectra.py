import numpy as np
from lib_QMout import *

eq_values = np.loadtxt("./energiesq.dat")
file_path = "./QM.out"
data = read_qmout(file_path)
total_stts, stts = getNstates(data)

multiplicidades = [i+1 for i, n in enumerate(stts) for _ in range(int(n))]
energias_SO = [E for E, m in zip(eq_values, multiplicidades) for _ in range(m)]
e_mat = getE(data, total_stts)
tmx,tmy,tmz = read_tm(data, total_stts)
dyson = read_Dyson(data, total_stts)

write_energy_differences(eq_values, prefix='energy_diff')
dyson = np.array(dyson)
write_SOC_diffs_withInt(energias_SO, dyson, prefix='soc_diff')

e_mat= np.array(e_mat)
for i in range(len(e_mat)):
    e_mat[i][2*i] = energias_SO[i]
write_NoSOC_QMout(data,stts,e_mat,tmx,tmy,tmz,dyson,"./test_qm.out")
