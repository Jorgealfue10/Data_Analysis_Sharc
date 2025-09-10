import numpy as np
from itertools import combinations

def read_qmout(file_path):
    with open(file_path,'r') as file:
        data=file.readlines()
    return data

def getNstates(data):
    onlyonce=False
    for line in data:
        if "states" in line and not onlyonce:
            parts=line.split()
            stts=parts[1:]
            onlyonce=True
            print(parts)

    total_stts=0
    for i in range(len(stts)):
        total_stts+=(i+1)*int(stts[i])
    
    return total_stts,stts

def getE(data,stts):
    start_readE=False
    e_mat=[]
    nsts=""+str(stts)+" "+str(stts)

    for line in data:
        if "Hamiltonian Matrix" in line:
            start_readE=True
            continue

        if start_readE:
            if nsts in line:
                continue
            elif line.strip() == "":
                break
            else:
                parts=line.split()
                e_mat.append(list(map(float,parts)))

    return e_mat

def map_acoplamientos(matrix1, matrix2, mult1, mult2, total_stts, final_matrix):
    """
    Mapea los valores fuera de la diagonal de dos matrices de acoplamiento en la matriz final.
    
    :param matrix1: np.array, primera matriz original de acoplamiento (N1 x 2N1, compleja)
    :param matrix2: np.array, segunda matriz original de acoplamiento (N2 x 2N2, compleja)
    :param mult1: list, lista con el número de estados por cada multiplicidad
    :param mult2: list, lista con el número de estados por cada multiplicidad
    :param total_stts: list, lista con el número total de estados por cada multiplicidad
    :param final_matrix: np.array, matriz destino con energías en la diagonal (M x 2M, compleja, M=N1+N2)
    :return: np.array, matriz final con acoplamientos correctamente mapeados
    """
    # Convertir listas a arrays de NumPy si es necesario
    matrix1 = np.array(matrix1)
    matrix2 = np.array(matrix2)
    final_matrix = np.array(final_matrix)

    def convertir_a_compleja(matrix):
        return matrix[:, ::2] + 1j * matrix[:, 1::2]
    
    complex_matrix1 = convertir_a_compleja(matrix1)
    complex_matrix2 = convertir_a_compleja(matrix2)
    complex_finalmat = convertir_a_compleja(final_matrix)
    
    # Checking same length of states.
    if len(total_stts) != len(mult1):
        diff=len(total_stts)-len(mult1)
        for i in range(diff):
            mult1.append(0)
    elif len(total_stts) != len(mult2):
        diff=len(total_stts)-len(mult2)
        for i in range(diff):
            mult2.append(0)

    M = final_matrix.shape[0]  # Tamaño de la matriz final
    mapped_matrix = np.zeros((M, M), dtype=complex)
    
    # Mapping the values from the original matrices to the final matrix
    # The indexes are mapped taking into account the #proj per multiplicities
    sidx,sidx1,sidx2=0,0,0
    for i in range(len(total_stts)):
        if total_stts[i] == mult1[i]:
            for j in range(int(total_stts[i])*(i+1)):
                for k in range(int(total_stts[i])*(i+1)):
                    new_j,new_k=sidx+j,sidx+k
                    if new_j != new_k:
                        mapped_matrix[new_j,new_k]=complex_matrix1[sidx1+j,sidx1+k]
                    else:
                        mapped_matrix[new_j,new_k]=complex_finalmat[new_j,new_k]
            sidx1+=int(total_stts[i])*(i+1)
        elif total_stts[i] == mult2[i]:
            for j in range(int(total_stts[i])*(i+1)):
                for k in range(int(total_stts[i])*(i+1)):
                    new_j,new_k=sidx+j,sidx+k
                    if new_j != new_k:
                        mapped_matrix[new_j,new_k]=complex_matrix2[sidx2+j,sidx2+k]
                    else:
                        mapped_matrix[new_j,new_k]=complex_finalmat[new_j,new_k]
            sidx2+=int(total_stts[i])*(i+1)
        sidx+=int(total_stts[i])*(i+1)

    return mapped_matrix

def read_tm(data,stts):
    tmx,tmy,tmz=[],[],[]

    start_read=False
    tmcoord=0
    nts=""+str(stts)+" "+str(stts)

    for line in data:
        if "Dipole Moment Matrices" in line:
            start_read=True
            continue

        if start_read:
            if nts in line:
                tmcoord+=1
                continue
            if tmcoord==1:
                parts=line.split()
                tmx.append(map(float,parts))
            elif tmcoord==2:
                parts=line.split()
                tmy.append(map(float,parts))
            elif tmcoord==3:
                if ("! 11 Prop" in line or line.strip() == "" or "! 6 Overlap" in line):
                    break
                parts=line.split()
                tmz.append(map(float,parts))
                
    return tmx,tmy,tmz

def read_Dyson(data,stts):
    dyson=[]
    start_read=False
    nsts=""+str(stts)+" "+str(stts)

    for line in data:
        if "Property Matrix" in line:
            start_read=True
            continue

        if start_read:
            if nsts in line:
                continue
            elif line.strip() == "":
                break
            else:
                parts=line.split()
                dyson.append(list(map(float,parts)))

    return dyson

def getxyz(data,nat):
    strt_xyz=False
    xyz=[]

    for line in range(len(data)):
        if line > 1 and not strt_xyz:
            strt_xyz=True

        if strt_xyz:
            if "unit" not in data[line]:
                xyz.append(data[line].split()[1:])
            else:
                break

    return xyz

def write_output(PHPHM_stts_num,final_mat,dyson_mat,tmx,tmy,tmz,outfile):
    # Writing the new output
    #with open("PHPHM_SOC.out","w") as file:
    with open(outfile,"w") as file:

        file.write("! 0 Basic information\n")
        file.write("states ")
        for i in range(len(PHPHM_stts_num)):
            file.write("{} ".format(PHPHM_stts_num[i]))
        file.write("\n")
        file.write("nmstates ")
        file.write("{}\n".format(final_mat.shape[0]))
        file.write("natom 2\n")
        file.write("npc 0\n")
        file.write("charges 0 0\n")
        file.write("\n")
        file.write("! 1 Hamiltonian Matrix ({}x{}, complex)\n".format(final_mat.shape[0], final_mat.shape[0]))
        file.write("{} {}\n".format(final_mat.shape[0], final_mat.shape[0]))

        for i in range(final_mat.shape[0]):
            for j in range(final_mat.shape[1]):
                file.write("{:18.12f} {:18.12f} ".format(final_mat[i,j].real,final_mat[i,j].imag))
            file.write("\n")
        
        file.write("\n")
        file.write("! 2 Dipole Moment Matrices (3x{}x{}, complex)\n".format(final_mat.shape[0], final_mat.shape[0]))
        file.write("{} {}\n".format(final_mat.shape[0], final_mat.shape[0]))

        for d in dyson_mat[:]:
            file.write(" ".join(f"{dy:.12f}" for dy in (list(d)))+"\n")

        file.write("{} {}\n".format(final_mat.shape[0], final_mat.shape[0]))
        for y in tmy:
            y=list(y)
            y[:]=[0]*len(y)
            file.write(" ".join(f"{ty:.12f}" for ty in (list(y)))+"\n")
        
        file.write("{} {}\n".format(final_mat.shape[0], final_mat.shape[0]))
        for z in tmz:
            z=list(z)
            z[:]=[0]*len(z)
            file.write(" ".join(f"{tz:.12f}" for tz in (list(z)))+"\n")

        file.write("\n")
        file.write("! 11 Property Matrix ({}x{}, complex)\n".format(len(list(dyson_mat)), len(list(dyson_mat[0]))))
        file.write("{} {}\n".format(len(list(dyson_mat)), len(list(dyson_mat[0]))))

        for x in tmx:
            file.write(" ".join(f"{tx:.12f}" for tx in (list(x)))+"\n")

        file.write("\n")

        file.write("! 8 Runtime\n")
        file.write("6.90000000000000E+001\n")

def write_NoSOC_QMout(data,stts,e_mat,tmx,tmy,tmz,dyson,file_out):
    
    with open(file_out,'w') as file:

        file.write("! 0 Basic information\n")
        file.write("states ")
        for i in range(len(stts)):
            file.write("{} ".format(stts[i]))
        file.write("\n")
        file.write("nmstates ")
        file.write("{}\n".format(e_mat.shape[0]))
        file.write("natom 2\n")
        file.write("npc 0\n")
        file.write("charges 0 0\n")
        file.write("\n")
        file.write("! 1 Hamiltonian Matrix ("+str(len(e_mat))+"x"+str(len(e_mat))+", complex)\n")
        file.write(str(len(e_mat))+" "+str(len(e_mat))+"\n")

        for e in e_mat:
            file.write(" ".join(f"{en:.12f}" for en in (list(e)))+"\n")

        file.write("\n")

        file.write("! 2 Dipole Moment Matrices (3x"+str(len(tmx))+"x"+str(len(tmx))+", complex)\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for y in tmy:
            y=list(y)
            y[:]=[0]*len(y)
            file.write(" ".join(f"{ty:.12f}" for ty in (list(y)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for d in dyson:
            file.write(" ".join(f"{dy:.12f}" for dy in (list(d)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for z in tmz:
            z=list(z)
            z[:]=[0]*len(z)
            file.write(" ".join(f"{tz:.12f}" for tz in (list(z)))+"\n")

        file.write("\n")
        file.write("! 11 Property Matrix ("+str(len(dyson))+"x"+str(len(dyson))+", complex)\n")
        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for x in tmx:
            file.write(" ".join(f"{tx:.12f}" for tx in (list(x)))+"\n")

        file.write("\n")

        file.write("! 8 Runtime\n")
        file.write("6.90000000000000E+001")

def write_NoSOC_NoDyson_QMout(data,stts,e_mat,tmx,tmy,tmz,dyson,file_out):
    
    with open(file_out,'w') as file:

        file.write("! 0 Basic information\n")
        file.write("states ")
        for i in range(len(stts)):
            file.write("{} ".format(stts[i]))
        file.write("\n")
        file.write("nmstates ")
        file.write("{}\n".format(e_mat.shape[0]))
        file.write("natom 2\n")
        file.write("npc 0\n")
        file.write("charges 0 0\n")
        file.write("\n")
        file.write("! 1 Hamiltonian Matrix ("+str(len(e_mat))+"x"+str(len(e_mat))+", complex)\n")
        file.write(str(len(e_mat))+" "+str(len(e_mat))+"\n")

        for e in e_mat:
            file.write(" ".join(f"{en:.12f}" for en in (list(e)))+"\n")

        file.write("\n")

        file.write("! 2 Dipole Moment Matrices (3x"+str(len(tmx))+"x"+str(len(tmx))+", complex)\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for y in tmy:
            y=list(y)
            y[:]=[0]*len(y)
            file.write(" ".join(f"{ty:.12f}" for ty in (list(y)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for d in dyson:
            file.write(" ".join(f"{dy:.12f}" for dy in (list(d)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for z in tmz:
            z=list(z)
            z[:]=[0]*len(z)
            file.write(" ".join(f"{tz:.12f}" for tz in (list(z)))+"\n")

        file.write("\n")
        file.write("! 11 Property Matrix ("+str(len(dyson))+"x"+str(len(dyson))+", complex)\n")
        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for x in tmx:
            file.write(" ".join(f"{tx:.12f}" for tx in (list(x)))+"\n")

        file.write("\n")

        file.write("! 8 Runtime\n")
        file.write("6.90000000000000E+001")

def write_energy_differences(energies, prefix='diff'):
    labels = [f"e{i+1}" for i in range(len(energies))]     
    for (i1, e1), (i2, e2) in combinations(enumerate(energies), 2):
        delta_e = abs(e2 - e1)
        label1 = labels[i1]
        label2 = labels[i2]
        filename = f"{prefix}_{label1}_{label2}.dat"
        with open(filename, 'w') as f:
            f.write(f"{delta_e:.8f}  1.0\n")

def write_SOC_diffs_withInt(energies, dyson, prefix='soc_diff'):
    labels = [f"e{i+1}" for i in range(len(energies))]     
    for i in range(len(energies)):
        for j in range(i+1, len(energies)):
            delta_e = abs(energies[j] - energies[i])
            label1 = labels[i]
            label2 = labels[j]
            filename = f"{prefix}_{label1}_{label2}.dat"
            intensity = delta_e * abs(dyson[i][2*j])**2
            with open(filename, 'w') as f:
                f.write(f"{delta_e:.8f}  {intensity:.8f}\n")