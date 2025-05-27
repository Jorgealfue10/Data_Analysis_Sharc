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

def write_MRCI_input(geoms,frozen,closed,occ,stts):
    """
    Escribe inputs de MRCI para todas las geometrías.

    :param geoms: array (nat x nicond x 3) con coordenadas
    :param frozen: núm. de orbitales congelados
    :param closed: núm. de orbitales cerrados
    :param occ: núm. de orbitales ocupados
    :param stts: lista de estados [[electrones, estados], ...]
    """

    for i in range(geoms.shape[1]):
        dir_path = makedirs("./MRCI/",i)
        with open(dir_path + "input.inp", 'w') as f:
            f.write("***,MOLPRO input for Davidson Energies MRCI \n")
            f.write("memory,300 \n")
            f.write("nosym \n bohrt \n")

            f.write("geometry={ \n")
            for j in range(geoms.shape[0]):
                element = "P" if j == 0 else "H"
                x,y,z = geoms[j,i]
                f.write(f" {element} {x:.6f} {y:.6f} {z:.6f} \n")
            f.write("} \n \n")

            f.write("gprint,orbitals,civectors; \n")
            f.write("gthresh,thrprint=0.,printci=0.0000000500; \n \n")

            f.write("basis=AV5Z \n")

            f.write("{hf \n wf,16,1,2 \n } \n")

            f.write("{multi, \n")
            f.write("frozen," + str(frozen) + "\n closed, " + str(closed) + " \n occ, " + str(occ) + " \n")
            for i, states in enumerate(stts):
                if isinstance(states, list) and states[1] != 0:
                    f.write(f"wf,{states[0]},1,{i} \n state,{states[1]} \n")
            f.write("} \n \n")

            ePH_vars = []
            ePHQ_vars = []

            for k, states in enumerate(stts):
                if isinstance(states, list) and states[1] != 0:
                    f.write("{mrci \n")
                    f.write("core,1 \n")
                    f.write(f"wf,{states[0]},1,{k} \n state,{states[1]} \n")
                    f.write("maxiter,250,1000 \n")
                    f.write("} \n")
                    if states[1] < 1:
                        f.write(f"ePH{k}=energy \n")
                        f.write(f"ePHQ{k}=energd \n")
                    else:
                        f.write(f"ePH{k}=energy(1) \n")
                        f.write(f"ePHQ{k}=energd(1) \n")

                    ePH_vars.append(f"ePH{k}")
                    ePHQ_vars.append(f"ePHQ{k}")
            
            f.write("table, " + ", ".join(ePH_vars) + "\n")
            f.write(f"save, energies.dat \n")

            f.write("table, " + ", ".join(ePHQ_vars) + "\n")
            f.write(f"save, energiesQ.dat \n")

            

data=read_file("initconds")
geoms=get_geoms(data, 2, 50)
for i in range(51):
    print(geoms[:,i,:])

sttsPH = [[16,0],[15,5],[16,1],[15,1]]
sttsPH2 = [[16,0],5,1,1]
sttsPH3 = [[16,0],5,1,1]

write_MRCI_input(geoms,0,1,10,sttsPH)