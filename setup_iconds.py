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
        f.write(f"Ninit \t {niconds}\n")
        f.write(f"Natoms \t {nat}\n")
        f.write(f"Reor \t None\n")
        f.write(f"Temp \t 0.00000\n")
        f.write(f"Eref \t 0.00000\n")
        f.write(f"Eharm \t 0.00000\n")

        f.write(f"\n")
        
        eqposP = [0.00000000,0.00000000,0.00000000]
        eqposH1 = [0.00000000,0.00000000,req]
        
        for i in range(nat):
            element = "P" if i == 0 else "H"
            zval = 15.0 if i == 0 else 1.0
            mass = 30.97376199842 if i == 0 else 1.00784016
            zeros = [0.00000000,0.00000000,0.00000000]
            f.write(f"Atom \t {element} \t {zval} \t {eqposP} \t {mass} \t {zeros} \n")

        f.write(f"\n")
        f.write(f"\n")
        
        for icond in range(niconds):
            rval = r0 + icond * dr
            posP = [0.00000000,0.00000000,0.00000000]
            posH1 = [0.00000000,0.00000000,rval]

            f.write(f"Index \t {icond+1} \n")
            f.write(f"Atoms \n")
            for i in range(nat):
                element = "P" if i == 0 else "H"
                x,y,z = posP if i == 0 else posH1
                f.write(f" {element} {x:.8f} {y:.8f} {z:.8f} {mass} {zeros} \n")
            f.write(f"States \n")
            f.write(f"Ekin \t 0.00000 \t a.u. \n")
            f.write(f"Epot_harm \t 0.00000 \t a.u. \n")
            f.write(f"Epot \t 0.00000 \t a.u. \n")
            f.write(f"Etot_harm \t 0.00000 \t a.u. \n")
            f.write(f"Etot \t 0.00000 \t a.u. \n")
        
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

def write_MRCI_input(geoms,frozen,closed,occ,stts,basis):
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
        with open(dir_path / "input.inp", 'w') as f:
            f.write("***,MOLPRO input for Davidson Energies MRCI \n")
            f.write("memory,600 \n")
            f.write("nosym \nbohr \n")

            f.write("geometry={ \n")
            for j in range(geoms.shape[0]):
                element = "P" if j == 0 else "H"
                x,y,z = geoms[j,i]
                f.write(f" {element} {x:.6f} {y:.6f} {z:.6f} \n")
            f.write("} \n \n")

            f.write(f"basis={basis} \n")

            f.write("{hf \n wf,17,1,1 \n } \n")

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

def write_MRCC_input(geoms,stts,basis):
    """
    Generates MRCC input files for each geometry configuration.

    :param geoms: numpy array of shape (nat, niconds+1, 3) containing atomic coordinates for each geometry configuration.
    :param stts: list of states, where each state is represented by a list containing the state information required for the MRCC input.
    """

    for i in range(geoms.shape[1]):
        dir_path = makedirs("./MRCC/",i)
        with open(dir_path / "input.inp", 'w') as f:
            f.write("***,MOLPRO input for MRCC \n")
            f.write("memory,300 \n")
            f.write("nosym \nbohr \n")

            f.write("geometry={ \n")
            for j in range(geoms.shape[0]):
                element = "P" if j == 0 else "H"
                x,y,z = geoms[j,i]
                f.write(f" {element} {x:.6f} {y:.6f} {z:.6f} \n")
            f.write("} \n \n")

            f.write("gprint,orbitals,civectors; \n")
            f.write("gthresh,thrprint=0.,printci=0.0000000500; \n \n")

            f.write(f"basis={basis} \n")

            ePH_vars = []

            for k, states in enumerate(stts):
                if isinstance(states, list) and states[1] != 0:
                    f.write("{hf \n")
                    f.write(f"wf,{states[0]},1,{k} \n")
                    f.write("} \n \n")

                    f.write("mrcc,method=ccsdt,dir=mrccdir \n")
                    f.write(f"ePH{k}=energy \n")

                    ePH_vars.append(f"ePH{k}")

            f.write("table, " + ", ".join(ePH_vars) + "\n")
            f.write(f"save, energies.dat \n")

# data=read_file("initconds")
nicond=50
nat=2
geom1D_gen(nat,nicond,2.7000,1.5000,3.5000)

# geoms=get_geoms(data, nat, nicond)
# for i in range(nicond+1):
#     print(geoms[:,i,:])

# basis="AV5Z"
# write_MRCI_input(geoms,0,6,15,sttsPH2,basis)
# write_MRCC_input(geoms,sttsPH2,basis)
