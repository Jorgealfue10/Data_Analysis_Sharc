import numpy as np

def read_molpro_output(file_path):
    with open(file_path,'r') as file:
        data = file.readlines()
    return data


def extract_freqcord(data):
    freqcord=[]
    start=False
    blank_line_seen=False
    nat=0
    for line in data:
        # the check below was incorrectly written as
        # `if 'Nr' and 'Atom' and 'Charge' in line:` which in Python only
        # verifies the last condition due to short circuit evaluation.
        # The intention is to start reading when all keywords are present.
        if 'Nr' in line and 'Atom' in line and 'Charge' in line:
            start=True
            blank_line_seen=False
            continue
        if start:
            if line.strip() == "" and not blank_line_seen:
                blank_line_seen=True
                continue
            elif line.strip() == "" and blank_line_seen:
                break
            parts=line.split()
            if len(parts) == 6:
                atom=parts[1]
                x,y,z=map(float,parts[3:])
                freqcord.append((atom,x,y,z))
                nat+=1
    return freqcord,nat

def extract_geometry(data):
    geometry=[]
    start=False
    blank_line_seen=False

    for line in data:
        if "Current geometry" in line:
            start=True
            continue

        if start:
            if line.strip() == "" and not blank_line_seen:
                blank_line_seen=True
                continue
            elif line.strip() == "" and blank_line_seen:
                break
            parts=line.split()
            if len(parts) == 4:
                atom=parts[0]
                x,y,z=map(float,parts[1:])
                geometry.append((atom,x,y,z))
    return geometry

def is_number(char):
    try:
        float(char)
        return True
    except:
        return False

def extract_freqs(data):
    freqs=[]
    low_freq=[]

    displ=[]
    low_disp=[]

    freq_start=False
    lowf_freq_start=False

    first_vib_seen=False
    sec_vib_seen=False

    for line in data:
        if 'Projecting out rotations and translations' in line:
            freq_start=True
            continue
        if freq_start:
            if not first_vib_seen and line.strip() == "":
                continue
            elif first_vib_seen and not sec_vib_seen and line.strip() == "":
                continue
            elif first_vib_seen and sec_vib_seen and line.strip() == "":
                break

            if "Vibration" in line and not first_vib_seen:
                first_vib_seen=True
                continue
            elif "Vibration" in line and first_vib_seen and not sec_vib_seen:
                sec_vib_seen=True
                continue

            parts=line.split()
            if first_vib_seen and not sec_vib_seen:
                if len(parts) == 2 and not is_number(parts[0]):
                   continue
                elif len(parts) == 2 and is_number(parts[1]):
                   low_freq.append(float(parts[1]))
            elif first_vib_seen and sec_vib_seen:
                if len(parts) == 2 and not is_number(parts[0]):
                   continue
                elif len(parts) == 2 and is_number(parts[1]):
                   freqs.append(float(parts[1]))

    return low_freq,freqs

def extract_normal_modes(data,nfreq,nat):
    
    ncoord=3
    nm_freq=np.zeros((nat, nfreq, ncoord))

    nm_freq_start=False
    # number of 5-mode blocks needed to store all frequencies
    nblocks = (nfreq + 4) // 5
    print(nfreq,nblocks,nfreq%5)

    blank_line=[]
    for i in range(nblocks):
        blank_line.append(False)

    print(nblocks)
    coord=0
    at=0
    block=0

    for line in data:
        if 'Normal Modes' in line:
            nm_freq_start=True
            continue

        if nm_freq_start:
            if line.strip() == "" and not all(blank_line):
                blank_line[block]=True
                continue
            elif line.strip() == "" and all(blank_line):
                break
            parts=line.split()
            print(parts)
            if "X" in line or "Y" in line or "Z" in line and is_number(parts[1]):
                if block+1 == nblocks:
                    nm_freq[at,block*5:nfreq,coord]=list(map(float,parts[1:]))
                else:
                    nm_freq[at,block*5:block*5+5,coord]=list(map(float,parts[1:]))

                coord+=1
                if coord == 3:
                    coord=0
                    at+=1
                if at == nat:
                    coord=0
                    at=0
                    block+=1

    return nm_freq

def extract_nm_lowf(data,nlowf,nat):
    ncoord=3
    nm_lowf=np.zeros((nat, nlowf, ncoord))

    nm_lowf_start=False

    # each block contains up to 5 modes
    nblocks = (nlowf + 4) // 5
    blank_nm_lowf=False

    blank_line=[]
    for i in range(nblocks):
        blank_line.append(False)

    coord=0
    at=0

    block=0
    for line2 in data:
        if 'Normal Modes of low' in line2:
            nm_lowf_start=True
            continue

        if nm_lowf_start:
            if line2.strip() == "" and not all(blank_line[:]):
                blank_line[block]=True
                continue
            elif line2.strip() == "" and all(blank_line[:]):
                break

            parts=line2.split()
            if "X" in line2 or "Y" in line2 or "Z" in line2 and is_number(parts[1]):
                if block+1 == nblocks:
                    nm_lowf[at,block*5:nlowf,coord]=list(map(float,parts[1:]))
                else:
                    nm_lowf[at,block*5:block*5+5,coord]=list(map(float,parts[1:]))

                coord+=1
                if coord == 3:
                    coord=0
                    at+=1
                if at == nat:
                    coord=0
                    at=0
                    block+=1

    return nm_lowf

def write_molden(geometry,freqcord,low_freq,freqs,nm_lowf,nm_freq,nat,output_file):
    with open(output_file,'w') as file:
        file.write("[Molden Format]\n")
        file.write("[Atoms] Angs\n")
        for i,(atom,x,y,z) in enumerate(geometry):
            file.write(f"{atom} {i+1} {atom} {x:.6f} {y:.6f} {z:.6f}\n")

        file.write("[FREQ]\n")
        for f in low_freq:
            file.write(f"{f:.6f}\n")
        for f in freqs:
            file.write(f"{f:.6f}\n")

        file.write(f"[FR-COORD]\n")
        for i,(atom,x,y,z) in enumerate(freqcord):
            file.write(str(geometry[i][0])+f" {x:.6f} {y:.6f} {z:.6f}\n")

        file.write("[FR-NORM-COORD]\n")
        vib=0
        for nm in range(nm_lowf.shape[1]):
            file.write(f"Vibration           {vib+1}\n")
            for at in range(nat):
#               print(nm_lowf[at,nm])
                file.write(" ".join(f"{coord:.8f}" for coord in nm_lowf[at,nm,:]) + "\n")
            vib+=1

        for nm in range(nm_freq.shape[1]):
            file.write(f"Vibration           {vib+1}\n")
            for at in range(nat):
#               print(nm_freq[at,nm])
                file.write(" ".join(f"{coord:.8f}" for coord in nm_freq[at,nm,:]) + "\n")
            vib+=1
