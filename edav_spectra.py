import numpy as np
from lib_QMout import *
from itertools import combinations

eq_values = np.loadtxt("./energiesq.dat")
file_path = "./QM.out"
data = read_qmout(file_path)
total_stts, stts = getNstates(data)

def write_energy_differences(energies, prefix='diff'):
    for e1, e2 in combinations(energies, 2):
        delta_e = abs(e2 - e1)
        # Nombres seguros (evitar puntos decimales raros o negativos)
        e1_label = f"{e1:.4f}".replace('.', 'p')
        e2_label = f"{e2:.4f}".replace('.', 'p')
        filename = f"{prefix}_{e1_label}_{e2_label}.dat"
        with open(filename, 'w') as f:
            f.write(f"{delta_e:.8f}  1.0\n")

multiplicidades = [i+1 for i, n in enumerate(stts) for _ in range(int(n))]
energias_SO = [E for E, m in zip(eq_values, multiplicidades) for _ in range(m)]
e_mat = getE(data, total_stts)
tmx,tmy,tmz = read_tm(data, total_stts)
dyson = read_Dyson(data, total_stts)

write_energy_differences(energias_SO, prefix='energy_diff')

e_mat= np.array(e_mat)
for i in range(len(e_mat)):
    e_mat[i][2*i] = energias_SO[i]
write_NoSOC_QMout(data,stts,e_mat,tmx,tmy,tmz,dyson,"./test_qm.out")
