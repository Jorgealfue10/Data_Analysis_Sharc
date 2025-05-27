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

def write_MRCI_input(geoms,frozen,closed,occ):
    """
    Escribe un archivo input para MRCI a partir de las geometrías proporcionadas.

    :param geoms: np.array, geometrías de las moléculas (nat x niconds x 3)
    """
    for i in range(geoms.shape[1]):
        dir_path = makedirs("./MRCI/",i)
        with open(dir_path + "input.inp", 'w') as f:
            f.write("***,MOLPRO input for Davidson Energies MRCI \n")
            f.write("memory,300 \n")
            f.write("nosym \n bohrt \n")

            f.write("geometry={ \n")
            for j in range(geoms.shape[0]):
                if j == 0:
                    f.write(f" P {geoms[j,i,0]} {geoms[j,i,1]} {geoms[j,i,2]}\n")
                else:   
                    f.write(f" H {geoms[j,i,0]} {geoms[j,i,1]} {geoms[j,i,2]}\n")
            f.write("} \n \n")

            f.write("gprint,orbitals,civectors; \n")
            f.write("gthresh,thrprint=0.,printci=0.0000000500; \n \n")

            f.write("basis=AV5Z \n")

            f.write("{hf \n wf,16,1,2 \n } \n")

            f.write("{multi, \n")
            f.write("frozen," + str(frozen) + "\n closed, " + str(closed) + " \n occ, " + str(occ) + " \n")
            f.write("wf,16,1,2 \n state,1 \n")
            f.write("wf,15,1,1 \n state,5 \n")
            f.write("wf,15,1,3 \n state,1 \n")
            f.write("} \n \n")

            f.write("{mrci \n")
            f.write("core,1 \n")
            f.write("wf,16,1,2 \n state,1 \n")
            f.write("maxiter,250,1000 \n")
            f.write("} \n")
            f.write("ePH=energy \n")
            f.write("ePHQ=energd \n")

            f.write("\n")

            f.write("{mrci \n")
            f.write("core,1 \n")
            f.write("wf,15,1,1 \n state,2 \n")
            f.write("maxiter,250,1000 \n")
            f.write("} \n")
            f.write("ePHM1=energy(1) \n")
            f.write("ePHM1Q=energd(1) \n")

            f.write("\n")

            f.write("{mrci \n")
            f.write("core,1 \n")
            f.write("wf,15,1,3 \n state,1 \n")
            f.write("maxiter,250,1000 \n")
            f.write("} \n")
            f.write("ePHM2=energy(1) \n")
            f.write("ePHM2Q=energd(1) \n")
            
            f.write("\n")

            f.write("table, (ePHM1-ePH)*27.2114, (ePHM2-ePH)*27.2114 \n")
            f.write("save, " + dir_path + "Ediff.dat \n")

            f.write("table, (ePHM1Q-ePHQ)*27.2114, (ePHM2Q-ePHQ)*27.2114 \n")
            f.write("save, " + dir_path + "EdiffQ.dat \n")
            

data=read_file("initconds")
geoms=get_geoms(data, 2, 50)
for i in range(51):
    print(geoms[:,i,:])

stts = [0,5,1,1]

write_MRCI_input(geoms,0,1,10)