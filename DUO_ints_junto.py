import numpy as np
import os ; import re
from pathlib import Path
import math ; from collections import defaultdict
import argparse

#-------------------------------------------------------------------------------#
# EXTRAER AUTOFUNCIONES VIBRACIONALES
#-------------------------------------------------------------------------------#
def is_head_eigenvib(line:str) -> bool:
    """
    Detecta la linea que marca el comienzo de 
    una nueva autofunción de un nivel vibracional
    Formato: 1 0.0000 1 0 X 3Sigmam -> 
    -------> Indice Energia Estado NVib Nombre Estado
    """

    s = line.strip()
    if not s:
        return False
    
    parts = s.split()
    #Si no tiene más de 2 campos
    if len(parts) < 2:
        return False
    
    #POR SI ACASOS
    #Si el primer campo no es un entero
    if not re.fullmatch(r"[+-]?\d+", parts[0]):
        return False
    
    #Si el segundo campo no es un valor decimal
    if not re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)([Ee][+-]?\d+)?", parts[1]):
        return False
    
    return True

#-------------------------------------------------------------------------------#
#Parsing Vibrational Eigenfunctions
#-------------------------------------------------------------------------------#

def parse_duo_vib_einfun(fname:str, npoints: int, nvib: int):
    """
    Parsea el archivo de autofunciones de vibración de DUO
    """
    f = open(fname, 'r')
    lines = f.readlines()

    vibmat = np.zeros((npoints,nvib))
    evect = np.zeros((nvib))
    ncurrpoint = 0 ; ncurrvib = 0
    for i,line in enumerate(lines):
        if is_head_eigenvib(line):
            ncurrpoint = 0
            ncurrvib = int(line.split()[3])
            if ncurrvib >= nvib:
                break
            evect[ncurrvib] = float(line.split()[1])
        elif "End of contracted basis" in line:
            break
        else:
            vibmat[ncurrpoint,ncurrvib] = float(line)
            ncurrpoint += 1

    return vibmat

#-------------------------------------------------------------------------------#
# Reading Dyson Matrix
#-------------------------------------------------------------------------------#

def dyson_mat(fname:str):
    dysvals = np.loadtxt(fname, dtype=np.float64)
    rvals = dysvals[:,0]; dys_N_C1 = dysvals[:,1]
    dysmat_N_C1 = np.zeros((dysvals.shape[0],dysvals.shape[0]))

    for i in range(dysvals.shape[0]):
            dysmat_N_C1[i][i] = dys_N_C1[i]*0.393456

    return rvals,dysmat_N_C1

#-------------------------------------------------------------------------------#
# Calculating BraKet
#-------------------------------------------------------------------------------#

def bra_ket(vib_ini, vib_fin, dymat):
    intensity = (vib_fin.T @ (dymat @ vib_ini))
    return intensity
def bk_Coeff(vib_ini, vib_fin, dymat, coef_ini, coef_fin):
    if np.array(coef_ini).shape[0] == 0 or np.array(coef_fin).shape[0] == 0:
        return bra_ket(vib_ini, vib_fin, dymat)
    
    coeffs = np.dot(np.array(coef_fin).conj(), np.array(coef_ini))
    intensity = (vib_fin.T @ (dymat @ vib_ini))
    intensity = intensity * coeffs
    return intensity

#-------------------------------------------------------------------------------#
# Get normalization factor Q(T) rotational levels including Projs
#-------------------------------------------------------------------------------#

def rot_normf(T, Jener, ZPE, Estt, numvib:int, numJ:int, Projs:int, index):

    Eh_to_cm = 219474.6
    kb = 3.166811563e-6  # Eh/K

    # calcular Emin global coherente
    Emin = 1e99
    for i in range(numvib):
        for j in range(numJ):
            for k in range(Projs):
                Ei = (Jener[i,j,k] + ZPE)/Eh_to_cm + Estt
                if Ei < Emin:
                    Emin = Ei

    # construir Q(T)
    Q = np.zeros((numvib))
    for i in range(numvib):
        for j in range(numJ):
            for k in range(Projs):
                Ei = (Jener[i,j,k] + ZPE)/Eh_to_cm + Estt
                Ei_rel = Ei - Emin
                jval=index[i][j][k][1]
                Q[i] += (2*jval+1) * np.exp(-Ei_rel/(kb*T))

    return Q

#-------------------------------------------------------------------------------#
# Get normalization factor Q(T) vibrational levels including Projs
#-------------------------------------------------------------------------------#

def vib_normf(T, Jener, ZPE, Estt, numvib:int, numJ:int, Projs:int, index):

    Eh_to_cm = 219474.6
    kb = 3.166811563e-6  # Eh/K

    # calcular Emin global coherente
    Emin = 1e99
    for i in range(numvib):
        Ei = (Jener[i][0][1] + ZPE)/Eh_to_cm + Estt
        if Ei < Emin:
            Emin = Ei

    # construir Q(T)
    Q = 0.0
    for i in range(numvib):
        Ei = (Jener[i][0][1] + ZPE)/Eh_to_cm + Estt
        Ei_rel = Ei - Emin
        Q += np.exp(-Ei_rel/(kb*T))

    return Q

#-------------------------------------------------------------------------------#
#Get Evals by Vib, J, Omega, Sigma, Lambda 
#-------------------------------------------------------------------------------#

def getJvals_compInd(fname:str, nvib: int, nj: int, mult: int, nLambda:int):
    f = open(fname, 'r')
    lines = f.readlines()

    S=(mult-1)/2
    nSigma = mult
    parity = 2
    jn = 0
    for i,line in enumerate(lines):
        if i == 0:
            if nLambda == 0:
                nOmega = mult
            else:
                nOmega = (nLambda+1) * mult
            nLind = nLambda + 1 
            Jvals = np.zeros((nvib,nj,nOmega))
            index_list = np.zeros((nvib,nj,nOmega,7))
            key_dict = np.zeros((nvib,nj,nOmega,6))

        parts = line.split()
        vibn = int(parts[4])
        jn = int(round(float(parts[0])))

        #Getting indexes for J projections
        mval = float(parts[8])
        if mval % 1 != 0 and nLambda != 0:
            midx = int(round(mval+S+1))
        else:
            midx = int(round(mval+S))

        if nLambda == 0:
            nL = 0
        else:
            if float(parts[5]) == -nLambda:
                nL = 0
            else:
                nL = 1
        nS = int(round(float(parts[7])+S))

        pm = parts[9]
        if pm == '-':
            p = 1
        else:
            p = 0
        
        if vibn < nvib:
            # print(i,vibn,float(parts[0]),float(parts[8]),float(parts[7]),float(parts[5]),p)
            Jvals[vibn,jn,midx] = float(parts[2])
            index_list[vibn,jn,midx,:] = vibn,float(parts[0]),float(parts[8]),float(parts[7]),float(parts[5]),p,int(parts[1]) # [0] Vib, [1] J, [2] m, [3] S, [4] Lambda, [5] Parity, [6] Index
            key_dict[vibn,jn,midx,:] = float(parts[0]),float(parts[8]),float(parts[7]),float(parts[5]),p,int(parts[1]) # [1] J, [2] m, [3] S, [4] Lambda, [5] Parity, [6] Index

    return Jvals,nOmega,index_list,key_dict

#-------------------------------------------------------------------------------#
#Kronecker delta
#-------------------------------------------------------------------------------#

def kronecker_delta(i,j):
    if i == j:
        return 1.0
    else:
        return 0.0

#-------------------------------------------------------------------------------#
#Factorial
#-------------------------------------------------------------------------------#
def fact(n):
    f=1
    if n%1 == 0:
        n=int(n)
    else:
        print(n)
        return 0
    for i in range(n):
        f=f*(i+1)
    return f 

#-------------------------------------------------------------------------------#
#Wigner 3j explicit expression
#-------------------------------------------------------------------------------#
def W3_exp(J1,J2,J3,m1,m2,m3):
    K=max(0,J2-J3-m1,J1-J3+m2)
    if K%1 != 0:
        return 0.0
    else: 
        K=int(K)
    N=min(J1+J2-J3,J1-m1,J2+m2)
    if N%1 != 0:
        return 0.0
    else: 
        N=int(N)
    
    deltaK = kronecker_delta(m1+m2+m3,0)
    m1pf = (-1)**(J1-J2-m3)

    #First element
    J1MJ2mJ3 = fact(J1+J2-J3)
    J1mJ2MJ3 = fact(J1-J2+J3)
    mJ1MJ2MJ3 = fact(-J1+J2+J3)
    J1MJ2MJ3M1 = fact(J1+J2+J3+1)
    first_root = math.sqrt((J1MJ2mJ3*J1mJ2MJ3*mJ1MJ2MJ3)/J1MJ2MJ3M1)

    #Second element
    J1mO1 = fact(J1-m1)
    J1MO1 = fact(J1+m1)
    J2mO2 = fact(J2-m2)
    J2MO2 = fact(J2+m2)
    J3mO3 = fact(J3-m3)
    J3MO3 = fact(J3+m3)
    second_root = math.isqrt(J1mO1*J1MO1*J2mO2*J2MO2*J3mO3*J3MO3)

    #Third element
    third_root = 0.0
    for k in range(K,N+1):
        m1pk = (-1)**k
        kf = fact(k)
        J1MJ2mJ3mk = fact(J1+J2-J3-k)
        J1mO1mk = fact(J1-m1-k)
        J2MO2mk = fact(J2+m2-k)
        J3mJ1MO1Mk = fact(J3-J2+m1+k)
        J3mJ1mO2Mk = fact(J3-J1-m2+k)
        coeff = m1pk/(kf*J1MJ2mJ3mk*J1mO1mk*J2MO2mk*J3mJ1MO1Mk*J3mJ1mO2Mk)
        third_root += coeff

    return deltaK*m1pf*first_root*second_root*third_root

#-------------------------------------------------------------------------------#
#Read coeff values according to index_list
#-------------------------------------------------------------------------------#

def read_coeff(filename):
    coef_dict = defaultdict(list)
    with open(filename) as f:
        for line in f:
            parts = line.split()

            if (parts[0]).isnumeric() == False: 
                continue
            else:
                vib = int(parts[5])
                J = float(parts[1])
                Om = float(parts[9])
                S = float(parts[8])
                L = float(parts[6])
                parity = int(parts[2])
                indx = int(parts[0])
                coeff = float(parts[3])

                key = (J,Om,S,L,parity,indx)
                coef_dict[key].append(coeff)

    for key in coef_dict:
        # print(coef_dict[key])
        coef_dict[key] = np.array(coef_dict[key])

    return coef_dict

#-------------------------------------------------------------------------------#
#Make intensity matrix including Rotational coeffs
#-------------------------------------------------------------------------------#

def intT_sep(PHener,PHMener,dysmat,PHvibs,PHMvibs,mask,Tvib,Trot,
                ZPEPH,ZPEPHM,PHsttsE,PHMsttsE,DJ:int,
                numvibPH:int,numvibPHM:int,numJPH:int,numJPHM:int,
                numOmPH:int,numOmPHM:int,indexPH,indexPHM,
                keyPH,keyPHM,coefPH,coefPHM):

    """
    Make the intensity matrix
    """
    Eh_to_cm = 219474.6
    dysmat = dysmat[mask] ; dysmat = dysmat[:,mask]
    relInt = np.zeros((numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM))
    evals = np.zeros((numvibPH,numvibPHM,numJPH,numJPHM,numOmPH,numOmPHM))

    normPHvib = vib_normf(Tvib,PHener,ZPEPH,PHsttsE,numvibPH,numJPH,numOmPH,indexPH)
    normPHrot = rot_normf(Trot,PHener,ZPEPH,PHsttsE,numvibPH,numJPH,numOmPH,indexPH)

    djlist = []
    initdj = -DJ
    for i in range(int(initdj+DJ),int(DJ+DJ)+1,1):
        djval = i -DJ
        djlist.append(djval)
    
    print(djlist)

    for i in range(numvibPH):
        for j in range(numvibPHM):

            PHvibsvn = PHvibs[:,i] ; PHvibsvn = PHvibsvn[mask]
            PHMvibsvn = PHMvibs[:,j] ; PHMvibsvn = PHMvibsvn[mask]

            for k in range(numJPH):
                for l in range(numJPHM):
                    for m in range(numOmPH):

                        PHv = (PHener[i,0,0]+ZPEPH)/Eh_to_cm
                        evib = PHv - ((np.min(PHener[:,0,0]) + ZPEPH)/Eh_to_cm)
                        vval=indexPH[i,k,m,0]
                        degvi = (np.exp(-evib/(3.166811563*10**(-6)*Tvib)))/normPHvib
                        PHv += PHsttsE

                        PHsumOmega = (PHener[i,k,m]+ZPEPH)/Eh_to_cm
                        eboltz = PHsumOmega - ((np.min(PHener[i,:,:]) + ZPEPH)/Eh_to_cm)
                        Jival=indexPH[i,k,m,1]
                        degJi = ((2*Jival+1)*np.exp(-eboltz/(3.166811563*10**(-6)*Trot)))/normPHrot[i]
                        PHsumOmega += PHsttsE

                        for n in range(numOmPHM):

                            keyPH_i = tuple(np.array(keyPH[i,k,m,:])) ; keyPHM_j = tuple(np.array(keyPHM[j,l,n,:]))
                            bk_val = bk_Coeff(PHvibsvn,PHMvibsvn,dysmat,coefPH[keyPH_i],coefPHM[keyPHM_j])

                            PHMsumOmega = ((PHMener[j,l,n])+ZPEPHM)/Eh_to_cm
                            PHMsumOmega += PHMsttsE

                            diff = abs(PHsumOmega - PHMsumOmega)

                            Jfval=indexPHM[j,l,n,1]

                            djval = Jfval - Jival
                            rotcoeff = 0.0
                            for dj in djlist:
                                rotcoeff += W3_exp(Jfval,DJ,Jival,indexPHM[j,l,n,2],dj,indexPH[i,k,m,2])
                            relInt[i,j,k,l,m,n] = degvi*degJi*((rotcoeff)*abs(bk_val))**2
                            evals[i,j,k,l,m,n] = diff
    return evals,relInt

#-------------------------------------------------------------------------------#
#Write files
#-------------------------------------------------------------------------------#

def dump_spectrum_long(filename, evals, relInt, index_listPH, index_listPHM,
                    energy_unit="eV", tol_I=0.0):
    """
    Guarda (E, I) + índices cuánticos discretos en formato long.

    Columns:
    vPH vPHM JPH JPHM OmPH OmPHM  E(unit)  I  idxPH  idxPHM
    """

    # Comprueba shapes básicas
    if evals.shape != relInt.shape:
        raise ValueError(f"evals.shape {evals.shape} != relInt.shape {relInt.shape}")

    if evals.ndim != 6:
        raise ValueError(f"Esperaba evals 6D, pero evals.ndim={evals.ndim}")

    nvPH, nvPHM, nJPH, nJPHM, nOmPH, nOmPHM = evals.shape

    Eh_to_eV = 27.2114 ; Eh_to_cm = 219474.63

    with open(filename, "w") as f:
        f.write(f"# columns: vPH JPH OmPH vPHM JPHM OmPHM  E({energy_unit})  I  \n")

        maxval = np.max(relInt[:,:,:,:,:,:])
        for i in range(nvPH):
            for j in range(nvPHM):
                for k in range(nJPH):
                    for l in range(nJPHM):
                        for m in range(nOmPH):
                            for n in range(nOmPHM):

                                I = relInt[i,j,k,l,m,n]#/maxval
                                if tol_I > 0.0 and abs(I) <= tol_I:
                                    continue

                                E = evals[i,j,k,l,m,n]

                                idxPH  = index_listPH[i,k,m,:]
                                idxPHM = index_listPHM[j,l,n,:]

                                for idx in idxPH:  # idxPHM: 
                                    f.write(f"{idx:.2f}  ")
                                for idx in idxPHM:
                                    f.write(f"{idx:.2f}  ")

                                if energy_unit == "eV":
                                    E *= Eh_to_eV
                                elif energy_unit == "cm^-1":
                                    E *= Eh_to_cm

                                f.write(f"{E: .10f} {I: .10f} \n")

#-------------------------------------------------------------------------------#
#Main
#-------------------------------------------------------------------------------#
def main():
    parser = argparse.ArgumentParser(
        description="Calcula valores de intensidad de espectro de fotoionización entre niveles rovibracionales."
    )
    parser.add_argument("-Tr", "--Trot", type=float, default=300.0, help="Temperatura rotacional en K")
    parser.add_argument("-Tv", "--Tvib", type=float, default=300.0, help="Temperatura vibracional en K")
    parser.add_argument("-ZPE", type=float, nargs=2, help="ZPE vibracional Sistema Neutro y catiónico (cm^-1)")
    parser.add_argument("-Etot", type=float, nargs=2, help="Energía total Sistema Neutro y Catiónico (Eh)")
    parser.add_argument("-nVJtot", type=int, nargs=2, help="Número de vibracionales y J totales Sistema Neutro y Catiónico")
    parser.add_argument("-DJ", type=float, default=1.5, help="DJ abs value")

    parser.add_argument("-Nvib", "--numvib", type=int, nargs=2, help="Número de vibracionales sistema neutro y catiónico")
    parser.add_argument("-NJ", "--numJ", type=int, nargs=2,help="Número de J sistema neutroy y catiónico")
    parser.add_argument("-mult", type=int, nargs=2, help="Multiplicidad sistema Neutroy y catiónico")
    parser.add_argument("-Lval", type=int, nargs=2, help="Valores Lambda sistema Neutro y Catiónico")

    parser.add_argument("-maskR", type=float, nargs=2, help="Rango de distancias")
    parser.add_argument("-npts", type=int, help="Número de puntos")

    parser.add_argument("-p", "--path", type=str, help="Ruta del directorio global")
    parser.add_argument("-dsys", type=str, nargs=2, help="Sistema Neutro y Catiónico")

    parser.add_argument("-o", "--output", type=str, default=None, help="Ruta del archivo de salida")

    args = parser.parse_args()

    Trot = args.Trot ; Tvib = args.Tvib
    nVib_tot, nJ_tot = args.nVJtot
    ZPEn = args.ZPE[0] ; ZPEc = args.ZPE[1]
    Etotn = args.Etot[0] ; Etotc = args.Etot[1]
    DJval = args.DJ
    numvibn = args.numvib[0] ; numvibc = args.numvib[1]
    numJn = args.numJ[0] ; numJc = args.numJ[1]
    multn = args.mult[0] ; multc = args.mult[1]
    Lvaln = args.Lval[0] ; Lvalc = args.Lval[1]

    pathor = args.path
    dsysn = args.dsys[0] ; dsysc = args.dsys[1]

    maskRi,maskRf = args.maskR
    npnts = args.npts

    fileout = args.output

    vibsn = parse_duo_vib_einfun(pathor+dsysn+"vibeigenvect_vib.chk", npnts, nVib_tot)
    vibsc = parse_duo_vib_einfun(pathor+dsysc+"vibeigenvect_vib.chk", npnts, nVib_tot)

    Jenern, nOmn, indexlsn, keysn = getJvals_compInd(pathor+dsysn+"rovibronic_energies.dat", nVib_tot, nJ_tot, multn, Lvaln)
    Jenerc, nOmc, indexlsc, keysc = getJvals_compInd(pathor+dsysc+"rovibronic_energies.dat", nVib_tot, nJ_tot, multc, Lvalc)

    coefn = read_coeff(pathor+dsysn+"vibeigenvect_vectors.chk")
    coefc = read_coeff(pathor+dsysc+"vibeigenvect_vectors.chk")

    rvals_Comp,dysmat = dyson_mat(pathor+"Dipole_moment_functions.dat")

    mask = (rvals_Comp > maskRi) & (rvals_Comp < maskRf)

    evals,relInt = intT_sep(Jenern,Jenerc,dysmat,vibsn,vibsc,mask,Tvib,Trot,
                            ZPEn,ZPEc,Etotn,Etotc,DJval,
                            numvibn,numvibc,numJn,numJc,nOmn,nOmc,
                            indexlsn,indexlsc,keysn,keysc,coefn,coefc)

    dump_spectrum_long(pathor+"/DJ"+str(DJval)+"/Tv"+str(Tvib)+"_Tr"+str(Trot)+fileout+".dat",
                        evals,relInt,indexlsn,indexlsc,energy_unit="eV",tol_I=0.0)

if __name__ == "__main__":
    main()

