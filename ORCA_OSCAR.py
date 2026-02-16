import numpy as np
import sys
import os
from pathlib import Path

def read_file(file_path):
    """
    Lee un archivo y devuelve su contenido como una lista de líneas.
    
    :param file_path: str, ruta del archivo a leer
    :return: list, contenido del archivo como lista de líneas
    """
    with open(file_path, 'r') as file:
        data = file.readlines()
    return data

def geom1D_gen(nat,niconds,req,r0,rf):
    """
    Genera un array para almacenar geometrías de una diatómica en formato ICONDs
    """

    dr=(rf - r0) / (niconds - 1)

    with open("initconds","w") as f:
        f.write(f"SHARC Initial Conditions File, version Jorge 2025\n")
        f.write(f"Ninit {niconds}\n")
        f.write(f"Natoms {nat}\n")
        f.write(f"Reor None\n")
        f.write(f"Temp 0.00000\n")
        f.write(f"Eref 0.00000\n")
        f.write(f"Eharm 0.00000\n")

        f.write(f"\n")
        
        eqposP = [0.00000000,0.00000000,0.00000000]
        eqposH1 = [0.00000000,0.00000000,req]
        
        f.write(f"Equilibrium \n")
        for i in range(nat):
            element = "P" if i == 0 else "H"
            zval = 15.0 if i == 0 else 1.0
            mass = 30.97376199842 if i == 0 else 1.00784016
            zeros = 0.00000000
            if i == 0:
                x,y,z = eqposP
                f.write(f"{element} {zval} {x:.6f} {y:.6f} {z:.6f} {mass} {zeros:6f} {zeros:6f} {zeros:6f} \n")
            else:
                x,y,z = eqposH1
                f.write(f"{element} {zval} {x:.6f} {y:.6f} {z:.6f} {mass} {zeros:6f} {zeros:6f} {zeros:6f} \n")

        f.write(f"\n")
        f.write(f"\n")
        
        for icond in range(niconds):
            rval = r0 + icond * dr
            posP = [0.00000000,0.00000000,0.00000000]
            posH1 = [0.00000000,0.00000000,rval]

            f.write(f"Index {icond+1} \n")
            f.write(f"Atoms \n")
            for i in range(nat):
                element = "P" if i == 0 else "H"
                mass = 30.97376199842 if i == 0 else 1.00784016
                x,y,z = posP if i == 0 else posH1
                f.write(f" {element} {x:.8f} {y:.8f} {z:.8f} {mass} {zeros:6f} {zeros:6f} {zeros:6f} \n")
            f.write(f"States \n")
            f.write(f"Ekin 0.00000 a.u. \n")
            f.write(f"Epot_harm 0.00000 a.u. \n")
            f.write(f"Epot 0.00000 a.u. \n")
            f.write(f"Etot_harm 0.00000 a.u. \n")
            f.write(f"Etot 0.00000 a.u. \n")
        
            f.write(f"\n")
            f.write(f"\n")


def get_geoms(data,nat,niconds):
    """
    Lee las geometrías en un archivo de QM.out y las devuelve en un array 3D
    
    :param data: list, contenido del archivo QM.out
    :param nat: int, número de átomos
    :param niconds: int, número de condiciones iniciales
    :return: array 3D, geometrías de las condiciones iniciales y la equilibrio
    """
    
    geoms=np.zeros((nat,niconds+1,3))
    for i,line in enumerate(data):
        print(line,i)
        if "Equilibrium" in line:
            curricond=0
            # print(line[i+1:i+2+nat],i)
            for j in range(nat):
                print(data[i+1+j],i+1+j)
                parts=data[i+1+j].split()
                x,y,z=float(parts[2]),float(parts[3]),float(parts[4])
                geoms[j,curricond,:]=[x,y,z]
            curricond+=1
        elif "Index" in line:
            for j in range(nat):
                parts=data[i+2+j].split()
                x,y,z=float(parts[2]),float(parts[3]),float(parts[4])
                geoms[j,curricond,:]=[x,y,z]
            curricond+=1
    return geoms

def makedirs(path,icond):
    """
    Crea un directorio para una condición inicial específica si no existe.
    
    :param path: str, ruta del directorio base
    :param icond: int, número de la condición inicial
    """

    dir_path = Path(path + f"/ICOND_{icond:05d}")
    dir_path.mkdir(parents=True, exist_ok=True)
    print(f"Directory created: {dir_path}" if os.path.exists(dir_path) else f"Directory already exists: {dir_path}")
    return dir_path

def orca_MRCI_input(geoms,nel,norb,stts,basis,chrge,multPC):
    """
    Escribe inputs de MRCI para todas las geometrías.

    :param geoms: array (nat x nicond x 3) con coordenadas
    :param frozen: núm. de orbitales congelados
    :param closed: núm. de orbitales cerrados
    :param occ: núm. de orbitales ocupados
    :param stts: lista de estados [[electrones, estados], ...]
    """

    for i in range(geoms.shape[1]):
        dir_path = makedirs("./",i)
        with open(dir_path / "input.inp", 'w') as f:
            f.write("!HF " + basis + " NoSym \n \n")

            f.write("%casscf \n")
            f.write("\t nel " + str(nel) + " \n")
            f.write("\t norb " + str(norb) + " \n")
            f.write("end \n \n")
            
            f.write("%mrci \n")
            f.write("\t CIType MRCI \n")
            f.write("\t DavidsonOpt Davidson1 \n")
            for i,nroots in enumerate(stts):
                if nroots != 0:
                    f.write("\t NewBlock "+str(i+1)+" 0\n")
                    f.write("\t \t NRoots "+str(nroots) + "\n")
                    if i%2 == 0:
                        f.write("\t \t Refs CAS(" + str(nel-1) + "," + str(norb) + ") end \n")
                    else:
                        f.write("\t \t Refs CAS(" + str(nel) + "," + str(norb) + ") end \n")
                    f.write("\t end \n ")
            f.write("end \n \n")
                
            f.write("* xyz "+ str(chrge) + " " + str(multPC) + " \n")
            for j in range(geoms.shape[0]):
                print(i)
                element = "P" if j == 0 else "C"
                x,y,z = geoms[j,i]
                f.write(f" {element} {x:.6f} {y:.6f} {z:.6f} \n")
            f.write("* \n \n")

data=read_file("initconds")
nicond=50
nat=2
geoms=get_geoms(data,nat,nicond)
nel=9 ; norb = 8
sttsOrcPC=[3,1,3,0,1]
basis="AV5Z"
charge=0
multPC=2
orca_MRCI_input(geoms,nel,norb,sttsOrcPC,basis,charge,multPC)
# write_MRCC_input(geoms,sttsPH2,basis)
