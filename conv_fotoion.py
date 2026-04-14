import numpy as np
import argparse

# --------------------------------------------------
# PERFILES
# --------------------------------------------------

def gaussian_profile(E, E0, sigma, I):
    return I * np.exp(-(E - E0)**2 / (2 * sigma**2))


def lorentzian_profile(E, E0, gamma, I):
    # Normalizada en área
    return I * (gamma / np.pi) / ((E - E0)**2 + gamma**2)


# --------------------------------------------------
# ESPECTRO GENERAL
# --------------------------------------------------

def build_spectrum(evals, intensities, E_grid, width, kind="gaussian", normalize=False):

    evals = np.asarray(evals).flatten()
    intensities = np.asarray(intensities).flatten()

    spectrum = np.zeros_like(E_grid)

    if kind == "gaussian":
        for E0, I0 in zip(evals, intensities):
            spectrum += gaussian_profile(E_grid, E0, width, I0)

    elif kind == "lorentzian":
        for E0, I0 in zip(evals, intensities):
            spectrum += lorentzian_profile(E_grid, E0, width, I0)

    else:
        raise ValueError("kind must be 'gaussian' or 'lorentzian'")

    if normalize:
        spectrum /= np.max(spectrum)

    return spectrum


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Spectrum convolution (clean version)")

    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)

    parser.add_argument("--emin", type=float, default=None)
    parser.add_argument("--emax", type=float, default=None)
    parser.add_argument("--dE", type=float, default=0.005)

    parser.add_argument("--width", type=float, default=0.05,
                        help="sigma (Gaussian) or gamma (Lorentzian) in eV")

    parser.add_argument("--kind", type=str, default="gaussian",
                        choices=["gaussian", "lorentzian"])

    parser.add_argument("-N", "--normalize", action="store_true")

    args = parser.parse_args()

    # ----------------------------------------------
    # Leer datos
    # ----------------------------------------------
    data = np.loadtxt(args.input, skiprows=1)

    # últimas dos columnas = energía + intensidad
    evals = data[:, -2]
    intensities = data[:, -1]

    # ----------------------------------------------
    # Grid de energía
    # ----------------------------------------------
    emin = args.emin if args.emin is not None else np.min(evals)
    emax = args.emax if args.emax is not None else np.max(evals)

    E_grid = np.arange(emin, emax, args.dE)

    # ----------------------------------------------
    # Construir espectro
    # ----------------------------------------------
    spectrum = build_spectrum(
        evals,
        intensities,
        E_grid,
        width=args.width,
        kind=args.kind,
        normalize=args.normalize
    )

    # ----------------------------------------------
    # Guardar
    # ----------------------------------------------
    np.savetxt(args.output, np.column_stack((E_grid, spectrum)),
            header="Energy(eV) Intensity")


if __name__ == "__main__":
    main()