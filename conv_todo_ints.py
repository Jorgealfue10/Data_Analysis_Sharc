import numpy as np
import os
import re
import glob
from collections import defaultdict


# =========================================================
# LORENTZIANA SOBRE GRILLA
# =========================================================
def lorentzian_spectrum_on_grid(evals, relInt, E_grid, gamma):

    energies = np.asarray(evals).flatten()
    intensities = np.asarray(relInt).flatten()

    spectrum = np.zeros_like(E_grid, dtype=float)

    for Ek, Ik in zip(energies, intensities):
        spectrum += Ik * gamma**2 / ((E_grid - Ek)**2 + gamma**2)

    return spectrum


# =========================================================
# AGRUPAR FICHEROS POR (Tv, Tr)
# =========================================================
def group_files_by_temp(folder):

    files = glob.glob(os.path.join(folder, "*.dat"))
    groups = defaultdict(list)

    for f in files:

        m = re.search(r"Tv([0-9.]+)_Tr([0-9.]+)", f)
        if not m:
            continue

        Tv = float(m.group(1))
        Tr = float(m.group(2))

        groups[(Tv, Tr)].append(f)

    return groups


# =========================================================
# LEER Y CONCATENAR DATOS
# =========================================================
def read_and_merge(files):

    all_data = []

    for f in files:
        data = np.loadtxt(f, skiprows=1)
        all_data.append(data)

    return np.vstack(all_data)


# =========================================================
# CALCULAR ESPECTROS
# =========================================================
def compute_spectra(data, E_grid, gamma):

    vN,JN,_,_,_,_,_,vC,JC,_,_,_,_,_,evals,relInt = data.T

    dV = (vC - vN)
    dJ = (JC - JN)

    spectra_dVdJ = {}

    # espectro total
    spec_total = lorentzian_spectrum_on_grid(evals, relInt, E_grid, gamma)

    # espectros parciales
    # for dv, dj in sorted(set(zip(dV, dJ))):

    #     mask = (dV == dv) & (dJ == dj)

    #     if np.sum(mask) == 0:
    #         continue

    #     spectra_dVdJ[(dv, dj)] = lorentzian_spectrum_on_grid(
    #         evals[mask],
    #         relInt[mask],
    #         E_grid,
    #         gamma
    #     )

    return spec_total #, spectra_dVdJ


# =========================================================
# GUARDAR RESULTADO
# =========================================================
def save_spectrum(outfile, E_grid, spec_total)#, spectra_dVdJ):

    cols = [E_grid]
    header = ["E"]

    # for dv, dj in sorted(spectra_dVdJ):
    #     cols.append(spectra_dVdJ[(dv, dj)])
    #     header.append(f"dV={dv}_dJ={dj}")

    # total al final
    cols.append(spec_total)
    header.append("TOTAL")

    data = np.column_stack(cols)

    np.savetxt(outfile, data, header=" ".join(header))


# =========================================================
# MAIN
# =========================================================
def main():

    folder = "DJ1.5"   # carpeta con los .dat
    dE = 0.001         # resolución en energía
    gamma = 0.01       # ancho lorentziano

    groups = group_files_by_temp(folder)

    print(f"Found {len(groups)} temperature groups\n")

    for (Tv, Tr), files in groups.items():

        print(f"Processing Tv={Tv}, Tr={Tr} ({len(files)} files)")

        # unir datos
        data = read_and_merge(files)

        evals = data[:,-2]

        emin = evals.min()
        emax = evals.max()

        E_grid = np.arange(emin, emax, dE)

        # calcular espectros
        spec_total, spectra_dVdJ = compute_spectra(data, E_grid, gamma)

        # nombre salida
        outfile = f"spectrum_Tv{Tv}_Tr{Tr}.dat"

        # guardar
        save_spectrum(outfile, E_grid, spec_total)#, spectra_dVdJ)

        print(f" -> saved {outfile}\n")


if __name__ == "__main__":
    main()