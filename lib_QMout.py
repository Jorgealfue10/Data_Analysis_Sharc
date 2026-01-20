import numpy as np
from itertools import combinations
from decimal import Decimal

# Reads QM.out file and returns its content as a list of lines
def read_qmout(file_path):
    with open(file_path,'r') as file:
        data=file.readlines()
    return data

# Gets the number of states and their multiplicities from the QM.out data
def getNstates(data):
    onlyonce=False
    for line in data:
        if "states" in line and not onlyonce:
            parts=line.split()
            stts=parts[1:]
            onlyonce=True
            # print(parts)

    total_stts=0
    for i in range(len(stts)):
        total_stts+=(i+1)*int(stts[i])
    
    return total_stts,stts

# Extracts the Hamiltonian matrix (energies) from the QM.out data
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

# Maps coupling values from two original matrices into the final matrix v1.0
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

# Reads transition moment matrices from QM.out data
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

# Reads Dyson matrix from QM.out data
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

# Gets the atomic coordinates from QM.out data
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

# Writes the final output QM.out file with SOC and changing the Dyson matrix
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
                cval_real = Decimal(final_mat[i,j].real)
                cval_imag = Decimal(final_mat[i,j].imag)
                print(cval_real,cval_imag)
                file.write("{:6.5e} {:6.5e} ".format(cval_real,cval_imag))
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

# Writes the final output QM.out file without SOC and changing Dyson matrix
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

# Writes the final output QM.out file without SOC and without Dyson matrix
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
        for x in tmx:
            x=list(x)
            file.write(" ".join(f"{tx:.12f}" for tx in (list(x)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for y in tmy:
            y=list(y)
            file.write(" ".join(f"{ty:.12f}" for ty in (list(y)))+"\n")

        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for z in tmz:
            z=list(z)
            file.write(" ".join(f"{tz:.12f}" for tz in (list(z)))+"\n")

        file.write("\n")
        file.write("! 11 Property Matrix ("+str(len(dyson))+"x"+str(len(dyson))+", complex)\n")
        file.write(str(len(tmx))+" "+str(len(tmx))+"\n")
        for d in dyson:
            file.write(" ".join(f"{dy:.12f}" for dy in (list(d)))+"\n")
        
        file.write("\n")

        file.write("! 8 Runtime\n")
        file.write("6.90000000000000E+001")

# Writes energy differences between states to separate files
def write_energy_differences(energies, prefix='diff'):
    labels = [f"e{i+1}" for i in range(len(energies))]     
    for (i1, e1), (i2, e2) in combinations(enumerate(energies), 2):
        delta_e = abs(e2 - e1)
        label1 = labels[i1]
        label2 = labels[i2]
        filename = f"{prefix}_{label1}_{label2}.dat"
        with open(filename, 'w') as f:
            f.write(f"{delta_e:.8f}  1.0\n")

# Writes energy differences and intensities based on Dyson matrix to separate files
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

# Transforms SHARC-type printed complex matrices to standard complex NumPy arrays and vice versa
def _to_complex(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M)
    assert M.ndim == 2 and M.shape[1] == 2*M.shape[0], f"Esperaba (N,2N), obtuve {M.shape}"
    return M[:, 0::2] + 1j*M[:, 1::2]

def _from_complex(C: np.ndarray) -> np.ndarray:
    C = np.asarray(C)
    N = C.shape[0]
    out = np.empty((N, 2*N), dtype=float)
    out[:, 0::2] = C.real
    out[:, 1::2] = C.imag
    return out

# Cumulative offsets for block indexing
def _cum_offsets(lengths: list[int]) -> list[int]:
    offs = [0]
    for x in lengths[:-1]:
        offs.append(offs[-1] + x)
    return offs

# Generalized mapping of coupling matrices into a final matrix, preserving diagonals
def map_acoplamientos_general(
    matrix1: np.ndarray,       # (N1 x 2N1) Re/Im intercalado
    matrix2: np.ndarray,       # (N2 x 2N2) Re/Im intercalado
    mult1: list[int],          # nº de ESTADOS por multiplicidad en matrix1 (sin MS)
    mult2: list[int],          # nº de ESTADOS por multiplicidad en matrix2 (sin MS)
    total_stts: list[int],     # nº TOTAL de ESTADOS por multiplicidad en la final (sin MS)
    final_matrix: np.ndarray   # (M x 2M) con energías en la diagonal (Re/Im intercalado)
) -> np.ndarray:
    """
    Mapea TODOS los acoplamientos (no diagonales) que ya están dentro de matrix1 y matrix2
    a la matriz final común, respetando multiplicidades y proyecciones MS, y preservando
    la diagonal (energías) de final_matrix.
    """
    L = len(total_stts)  # nº de multiplicidades (S=0, 1/2, 1, ...)
    # Alinear longitudes por si mult1/mult2 son más cortas:
    mult1 = list(mult1) + [0]*max(0, L - len(mult1))
    mult2 = list(mult2) + [0]*max(0, L - len(mult2))
    print(total_stts, mult1, mult2)
    # Proyecciones por multiplicidad: P = (m_idx+1) * nº_estados
    Ptot = [(m+1)*total_stts[m] for m in range(L)]
    P1   = [(m+1)*mult1[m]      for m in range(L)]
    P2   = [(m+1)*mult2[m]      for m in range(L)]
    print(Ptot, P1, P2)
    # Convertir a complejo (si la matriz existe)
    C1 = _to_complex(matrix1) if matrix1 is not None and np.size(matrix1) else None
    C2 = _to_complex(matrix2) if matrix2 is not None and np.size(matrix2) else None
    Cf = _to_complex(final_matrix)  # base: contiene la diagonal (energías)

    # Checks de tamaños globales
    Nf = Cf.shape[0]
    print(Nf,Ptot)
    assert sum(Ptot) == Nf, f"Suma proyecciones finales {sum(Ptot)} != Nf {Nf}"
    if C1 is not None:
        assert sum(P1) == C1.shape[0], f"Suma P1 {sum(P1)} != N1 {C1.shape[0]}"
    if C2 is not None:
        assert sum(P2) == C2.shape[0], f"Suma P2 {sum(P2)} != N2 {C2.shape[0]}"

    # Offsets por multiplicidad en la final:
    off_fin = _cum_offsets(Ptot)
    off1_fin = off_fin[:]                            # matrix1 empieza aquí dentro de cada mult
    off2_fin = [off_fin[m] + P1[m] for m in range(L)]# matrix2 va detrás dentro de cada mult

    # Offsets por multiplicidad en cada origen:
    off1_src = _cum_offsets(P1)
    off2_src = _cum_offsets(P2)

    # Partimos de Cf para preservar energías
    dest = Cf.copy()

    def copy_all_blocks(Cs, P, off_src, off_dst):
        """
        Copia TODOS los subbloques de Cs (incluyendo entre multiplicidades)
        a su posición en 'dest'. En bloques m==n, restaura su diagonal local.
        """
        if Cs is None:
            return
        for m in range(L):
            rm = P[m]
            if rm == 0: 
                continue
            rs = slice(off_src[m], off_src[m] + rm)
            rd = slice(off_dst[m], off_dst[m] + rm)
            for n in range(L):
                cn = P[n]
                if cn == 0:
                    continue
                cs = slice(off_src[n], off_src[n] + cn)
                cd = slice(off_dst[n], off_dst[n] + cn)
                block = Cs[rs, cs]
                if m == n:
                    keep = dest[rd, cd].copy()
                    dest[rd, cd] = block
                    k = min(rm, cn)
                    idx = np.arange(k)
                    dest[rd, cd][idx, idx] = keep[idx, idx]  # preserva la diagonal local
                else:
                    dest[rd, cd] = block

    # Copiamos todo lo que haya en matrix1 y matrix2
    copy_all_blocks(C1, P1, off1_src, off1_fin)
    copy_all_blocks(C2, P2, off2_src, off2_fin)

    # Asegurar la diagonal global (por si acaso)
    np.fill_diagonal(dest, np.diag(Cf))

    # Volvemos a (M x 2M) Re/Im intercalado
    return _from_complex(dest)
