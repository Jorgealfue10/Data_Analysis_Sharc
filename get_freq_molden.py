import numpy as np

molpro_out="MOLPRO.out"
molden_file="freq.molden"

data=read_molpro_output(molpro_out)
freqcord,nat=extract_freqcord(data)
geometry=extract_geometry(data)
low_freq,freqs=extract_freqs(data)
print(len(low_freq),len(freqs))
nm_freq=extract_normal_modes(data,len(freqs),nat)
nm_lowf=extract_nm_lowf(data,len(low_freq),nat)
write_molden(geometry,freqcord,low_freq,freqs,nm_lowf,nm_freq,nat,molden_file)
