import numpy as np
import argparse

# --------------------------------------------------
# PERFILES NORMALIZADOS EN ÁREA
# --------------------------------------------------

def gaussian_profile(E, E0, sigma, I):
    """
    Perfil gaussiano normalizado en área:
    integral = I
    """
    pref = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    return I * pref * np.exp(-(E - E0)**2 / (2.0 * sigma**2))


def lorentzian_profile(E, E0, gamma, I):
    """
    Perfil lorentziano normalizado en área:
    integral = I
    """
    return I * (gamma / np.pi) / ((E - E0)**2 + gamma**2)


# --------------------------------------------------
# ESPECTRO GENERAL
# --------------------------------------------------

def build_spectrum(
    evals,
    intensities,
    E_grid,
    width,
    kind="gaussian",
    normalize=False
):
    """
    Construye el espectro convolucionado.

    width:
        gaussian   -> sigma
        lorentzian -> gamma
    """

    evals = np.asarray(evals, dtype=float).flatten()
    intensities = np.asarray(intensities, dtype=float).flatten()

    spectrum = np.zeros_like(E_grid, dtype=float)

    # --------------------------------------------------
    # Gaussian
    # --------------------------------------------------
    if kind == "gaussian":

        sigma = width
        cutoff = 5.0 * sigma   # suficiente para gaussiana

        for E0, I0 in zip(evals, intensities):

            mask = np.abs(E_grid - E0) <= cutoff

            if not np.any(mask):
                continue

            spectrum[mask] += gaussian_profile(
                E_grid[mask],
                E0,
                sigma,
                I0
            )

    # --------------------------------------------------
    # Lorentzian
    # --------------------------------------------------
    elif kind == "lorentzian":

        gamma = width
        cutoff = 20.0 * gamma  # cola más larga

        for E0, I0 in zip(evals, intensities):

            mask = np.abs(E_grid - E0) <= cutoff

            if not np.any(mask):
                continue

            spectrum[mask] += lorentzian_profile(
                E_grid[mask],
                E0,
                gamma,
                I0
            )

    else:
        raise ValueError("kind must be 'gaussian' or 'lorentzian'")

    # --------------------------------------------------
    # Normalización opcional
    # --------------------------------------------------
    if normalize and np.max(spectrum) > 0.0:
        spectrum /= np.max(spectrum)

    return spectrum


# --------------------------------------------------
# MAIN
# --------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description="Proper stick-spectrum convolution"
    )

    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True
    )

    parser.add_argument(
        "--emin",
        type=float,
        default=None
    )

    parser.add_argument(
        "--emax",
        type=float,
        default=None
    )

    parser.add_argument(
        "--dE",
        type=float,
        default=None,
        help="Energy step of final grid"
    )

    parser.add_argument(
        "--width",
        type=float,
        default=0.05,
        help="sigma (Gaussian) or gamma (Lorentzian) in eV"
    )

    parser.add_argument(
        "--kind",
        type=str,
        default="gaussian",
        choices=["gaussian", "lorentzian"]
    )

    parser.add_argument(
        "-N", "--normalize",
        action="store_true"
    )

    args = parser.parse_args()

    # --------------------------------------------------
    # Leer datos
    # --------------------------------------------------
    data = np.loadtxt(args.input, skiprows=1)

    # últimas dos columnas = energía + intensidad
    evals = data[:, -2]
    intensities = data[:, -1]

    # --------------------------------------------------
    # Width
    # --------------------------------------------------
    width = float(args.width)

    # --------------------------------------------------
    # dE inteligente
    #
    # regla segura:
    # al menos ~10 puntos por ancho
    # --------------------------------------------------
    if args.dE is None:
        dE = width / 10.0
    else:
        dE = args.dE

    # --------------------------------------------------
    # Límites de energía
    # ampliados para no cortar colas
    # --------------------------------------------------
    if args.kind == "gaussian":
        extra = 5.0 * width
    else:
        extra = 20.0 * width

    emin = args.emin if args.emin is not None else np.min(evals) - extra
    emax = args.emax if args.emax is not None else np.max(evals) + extra

    E_grid = np.arange(emin, emax + dE, dE)

    # --------------------------------------------------
    # Construcción
    # --------------------------------------------------
    spectrum = build_spectrum(
        evals=evals,
        intensities=intensities,
        E_grid=E_grid,
        width=width,
        kind=args.kind,
        normalize=args.normalize
    )

    # --------------------------------------------------
    # Guardar
    # --------------------------------------------------
    np.savetxt(
        args.output,
        np.column_stack((E_grid, spectrum)),
        header="Energy(eV) Intensity"
    )


if __name__ == "__main__":
    main()