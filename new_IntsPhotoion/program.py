import argparse
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d
from tqdm import tqdm


EH_TO_CM = 219474.6
EH_TO_EV = 27.2114
BOHR_TO_ANGSTROM = 0.529177210903


# =============================================================================
# DUO level object
# =============================================================================

@dataclass(frozen=True)
class DuoLevel:
    system: str
    v: int
    J: float
    Omega: float
    Sigma: float
    Lambda: float
    parity: int
    duo_index: int
    energy_cm: float


@dataclass(frozen=True)
class CoeffComponent:
    coeff: complex
    v_basis: int
    Lambda: float
    Sigma: float
    Omega: float


# =============================================================================
# Small utilities
# =============================================================================

def as_list(x):
    if isinstance(x, (list, tuple, np.ndarray)):
        return list(x)
    return [x]


def resolve_system_file(base_path, system, filename):
    """
    Accepts both directory-style and prefix-style inputs.

    Directory style:
        base/PHX3Sm/vibeigenvect_vib.chk

    Prefix style:
        base/PHX3Smvibeigenvect_vib.chk
    """

    base_path = Path(base_path)

    candidate_dir = base_path / system / filename
    if candidate_dir.exists():
        return candidate_dir

    candidate_prefix = base_path / f"{system}{filename}"
    if candidate_prefix.exists():
        return candidate_prefix

    raise FileNotFoundError(
        f"Could not find {filename} for system {system}. Tried:\n"
        f"  {candidate_dir}\n"
        f"  {candidate_prefix}"
    )


def build_dyson_index_to_omega():
    """
    Maps Dyson file indices to physical projections.

    These indices are NOT DUO rovibronic level indices. They are the indices
    appearing in filenames like dyson_12_02.dat.
    """

    return {
        "PHX3Sm": {
            12: -1.0,
            13:  0.0,
            14:  1.0,
        },

        # Singlet Delta: the same Dyson file is used for both ±2 components.
        "PHa1D": {
            1: [-2.0, 2.0],
        },

        "PHMX2P": {
            2: -1.5,
            3: -0.5,
            7:  0.5,
            8:  1.5,
        },

        "PHMa4Sm": {
            15: -1.5,
            16: -0.5,
            17:  0.5,
            18:  1.5,
        },

        "PHMA2D": {
            5:  -2.5,
            6:  -1.5,
            10:  1.5,
            11:  2.5,
        },

        "PHM12Sm": {
            4: -0.5,
            9:  0.5,
        },
    }


# =============================================================================
# Vibrational eigenfunctions
# =============================================================================

def is_head_eigenvib(line: str) -> bool:
    s = line.strip()
    if not s:
        return False

    parts = s.split()
    if len(parts) < 2:
        return False

    if not re.fullmatch(r"[+-]?\d+", parts[0]):
        return False

    if not re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)([Ee][+-]?\d+)?", parts[1]):
        return False

    return True


def parse_duo_vib_einfun(fname, npoints: int, nvib: int):
    """
    Parses DUO vibrational contracted eigenfunctions.
    Returns an array with shape (npoints, nvib).
    """

    fname = Path(fname)
    vibmat = np.zeros((npoints, nvib), dtype=float)

    with open(fname, "r") as f:
        lines = f.readlines()

    ncurrpoint = 0
    ncurrvib = None

    for line in lines:
        if is_head_eigenvib(line):
            ncurrpoint = 0
            ncurrvib = int(line.split()[3])
            if ncurrvib >= nvib:
                break

        elif "End of contracted basis" in line:
            break

        else:
            if ncurrvib is None:
                continue
            if ncurrpoint >= npoints:
                raise ValueError(
                    f"Too many points while reading vib={ncurrvib} in {fname}"
                )
            vibmat[ncurrpoint, ncurrvib] = float(line)
            ncurrpoint += 1

    return vibmat


# =============================================================================
# DUO rovibronic levels and coefficients
# =============================================================================

def read_duo_levels(fname, system, nvib_max=None, j_max=None):
    """
    Reads DUO rovibronic_energies.dat as a flat list of real levels.
    No dense zero-filled arrays are used.
    """

    fname = Path(fname)
    levels = []

    with open(fname, "r") as f:
        for iline, line in enumerate(f, start=1):
            parts = line.split()

            if len(parts) < 10:
                continue

            J = float(parts[0])
            duo_index = int(parts[1])
            energy_cm = float(parts[2])
            v = int(parts[4])
            Lambda = float(parts[5])
            Sigma = float(parts[7])
            Omega = float(parts[8])
            parity = 1 if parts[9] == "-" else 0

            if nvib_max is not None and v >= nvib_max:
                continue

            if j_max is not None and J >= j_max:
                continue

            if duo_index <= 0:
                raise ValueError(
                    f"Invalid DUO index in {fname}, line {iline}: {duo_index}"
                )

            levels.append(
                DuoLevel(
                    system=system,
                    v=v,
                    J=J,
                    Omega=Omega,
                    Sigma=Sigma,
                    Lambda=Lambda,
                    parity=parity,
                    duo_index=duo_index,
                    energy_cm=energy_cm,
                )
            )

    return levels


def level_state_key(level):
    return (
        level.J,
        level.Omega,
        level.Sigma,
        level.Lambda,
        level.parity,
        level.duo_index,
    )


def read_coefficients(filename):
    filename = Path(filename)
    tmp = defaultdict(list)

    with open(filename, "r") as f:
        for line in f:
            parts = line.split()

            if len(parts) < 10:
                continue

            try:
                duo_index = int(parts[0])
                J = float(parts[1])
                parity = int(parts[2])
                coeff = complex(float(parts[3]), 0.0)
                v_basis = int(parts[5])
                Lambda = float(parts[6])
                Sigma = float(parts[8])
                Omega = float(parts[9])
            except ValueError:
                continue

            key = (
                J,
                Omega,
                Sigma,
                Lambda,
                parity,
                duo_index,
            )

            tmp[key].append(
                CoeffComponent(
                    coeff=coeff,
                    v_basis=v_basis,
                    Lambda=Lambda,
                    Sigma=Sigma,
                    Omega=Omega,
                )
            )

    return dict(tmp)


# =============================================================================
# Overlaps
# =============================================================================

def bra_ket_weighted(vib_ini, vib_fin, dyson_values):
    """
    Computes <chi_f | D(R) | chi_i> without building diag(D).
    """

    return np.dot(vib_fin, dyson_values * vib_ini)




# =============================================================================
# Wigner 3j explicit expression
# =============================================================================

def kronecker_delta(i, j):
    return 1.0 if i == j else 0.0


def fact(n):
    if n % 1 != 0:
        return 0

    n = int(n)
    if n < 0:
        return 0

    f = 1
    for i in range(n):
        f *= i + 1
    return f


def W3_exp(J1, J2, J3, m1, m2, m3):
    K = max(0, J2 - J3 - m1, J1 - J3 + m2)
    if K % 1 != 0:
        return 0.0
    K = int(K)

    N = min(J1 + J2 - J3, J1 - m1, J2 + m2)
    if N % 1 != 0:
        return 0.0
    N = int(N)

    deltaK = kronecker_delta(m1 + m2 + m3, 0)
    if deltaK == 0.0:
        return 0.0

    m1pf = (-1) ** int(J1 - J2 - m3)

    J1MJ2mJ3 = fact(J1 + J2 - J3)
    J1mJ2MJ3 = fact(J1 - J2 + J3)
    mJ1MJ2MJ3 = fact(-J1 + J2 + J3)
    J1MJ2MJ3M1 = fact(J1 + J2 + J3 + 1)

    if J1MJ2MJ3M1 == 0:
        return 0.0

    first_root = math.sqrt(
        (J1MJ2mJ3 * J1mJ2MJ3 * mJ1MJ2MJ3) / J1MJ2MJ3M1
    )

    J1mO1 = fact(J1 - m1)
    J1MO1 = fact(J1 + m1)
    J2mO2 = fact(J2 - m2)
    J2MO2 = fact(J2 + m2)
    J3mO3 = fact(J3 - m3)
    J3MO3 = fact(J3 + m3)

    second_arg = J1mO1 * J1MO1 * J2mO2 * J2MO2 * J3mO3 * J3MO3
    if second_arg < 0:
        return 0.0

    second_root = math.sqrt(second_arg)

    third_sum = 0.0
    for k in range(K, N + 1):
        kf = fact(k)
        J1MJ2mJ3mk = fact(J1 + J2 - J3 - k)
        J1mO1mk = fact(J1 - m1 - k)
        J2MO2mk = fact(J2 + m2 - k)
        J3mJ1MO1Mk = fact(J3 - J2 + m1 + k)
        J3mJ1mO2Mk = fact(J3 - J1 - m2 + k)

        denom = (
            kf
            * J1MJ2mJ3mk
            * J1mO1mk
            * J2MO2mk
            * J3mJ1MO1Mk
            * J3mJ1mO2Mk
        )

        if denom == 0:
            continue

        third_sum += ((-1) ** k) / denom

    return deltaK * m1pf * first_root * second_root * third_sum


# =============================================================================
# Dyson reading and splines
# =============================================================================

def read_dyson_raw(d_neutral, d_cation, dyson_path, bohr_to_angstrom=True):
    dyson_path = Path(dyson_path)
    raw = {}
    r_ref = None

    for idx_i, omega_i_vals in d_neutral.items():
        for idx_f, omega_f_vals in d_cation.items():
            fname = dyson_path / f"dyson_{idx_i:02d}_{idx_f:02d}.dat"

            if not fname.exists():
                raise FileNotFoundError(f"Missing Dyson file: {fname}")

            r, dyson = np.loadtxt(fname, unpack=True)

            if r_ref is None:
                r_ref = r.copy()
            else:
                if len(r) != len(r_ref) or not np.allclose(r, r_ref):
                    raise ValueError(f"Inconsistent R grid in {fname}")

            r_use = r * BOHR_TO_ANGSTROM if bohr_to_angstrom else r.copy()

            raw[(idx_i, idx_f)] = {
                "idx_i": idx_i,
                "idx_f": idx_f,
                "omega_i_values": as_list(omega_i_vals),
                "omega_f_values": as_list(omega_f_vals),
                "r": r_use,
                "dyson": np.asarray(dyson, dtype=float),
                "file": str(fname),
            }

    return raw


def combine_dyson_by_omega(raw, mode="quadrature"):
    if mode not in ["quadrature", "sum"]:
        raise ValueError("mode must be 'quadrature' or 'sum'")

    accum = defaultdict(lambda: None)
    sources = defaultdict(list)
    r_by_key = {}

    for (idx_i, idx_f), item in raw.items():
        r = item["r"]
        dyson = item["dyson"]

        for omega_i in item["omega_i_values"]:
            for omega_f in item["omega_f_values"]:
                key = (omega_i, omega_f)

                if key not in r_by_key:
                    r_by_key[key] = r
                else:
                    if len(r) != len(r_by_key[key]) or not np.allclose(r, r_by_key[key]):
                        raise ValueError(f"Non-matching R grids for Dyson key {key}")

                sources[key].append(item["file"])

                if accum[key] is None:
                    accum[key] = dyson**2 if mode == "quadrature" else dyson.copy()
                else:
                    accum[key] += dyson**2 if mode == "quadrature" else dyson

    combined = {}
    for key, val in accum.items():
        dyson_eff = np.sqrt(val) if mode == "quadrature" else val
        combined[key] = {
            "omega_i": key[0],
            "omega_f": key[1],
            "r": r_by_key[key],
            "dyson": dyson_eff,
            "sources": sources[key],
        }

    return combined


def make_dyson_splines(combined, spline_type="cubic", extrapolate=False):
    if spline_type not in ["cubic", "linear"]:
        raise ValueError("spline_type must be 'cubic' or 'linear'")

    splines = {}

    for key, item in combined.items():
        r = np.asarray(item["r"], dtype=float)
        dyson = np.asarray(item["dyson"], dtype=float)

        order = np.argsort(r)
        r = r[order]
        dyson = dyson[order]

        if np.any(~np.isfinite(r)) or np.any(~np.isfinite(dyson)):
            raise ValueError(f"Non-finite values in Dyson data for key {key}")

        if np.any(np.diff(r) <= 0):
            raise ValueError(f"R grid is not strictly increasing for key {key}")

        if spline_type == "cubic":
            splines[key] = CubicSpline(
                r,
                dyson,
                bc_type="natural",
                extrapolate=extrapolate,
            )
        else:
            splines[key] = interp1d(
                r,
                dyson,
                kind="linear",
                bounds_error=not extrapolate,
                fill_value="extrapolate" if extrapolate else np.nan,
            )

    return splines


def read_dyson_splines(
    d_neutral,
    d_cation,
    dyson_path,
    combine_mode="quadrature",
    spline_type="cubic",
    bohr_to_angstrom=True,
    extrapolate=False,
):
    raw = read_dyson_raw(
        d_neutral=d_neutral,
        d_cation=d_cation,
        dyson_path=dyson_path,
        bohr_to_angstrom=bohr_to_angstrom,
    )

    combined = combine_dyson_by_omega(raw, mode=combine_mode)

    splines = make_dyson_splines(
        combined=combined,
        spline_type=spline_type,
        extrapolate=extrapolate,
    )

    return splines, raw, combined


# =============================================================================
# Transition table
# =============================================================================

def rotational_factor(J_i, J_f, Omega_i, Omega_f, K_values):
    total = 0.0

    for K in K_values:
        if not (abs(J_i - J_f) <= K <= (J_i + J_f)):
            continue

        if not np.isclose((J_i + J_f + K) % 1, 0.0):
            continue

        for dOmega in np.arange(-K, K + 1.0, 1.0):
            if not np.isclose(Omega_i + Omega_f - dOmega, 0.0):
                continue

            total += W3_exp(
                J_i,
                J_f,
                K,
                Omega_i,
                Omega_f,
                -dOmega,
            )

    return total


def transition_energy_Eh(level_i, level_f, ZPE_i_cm, ZPE_f_cm, Eelec_i_Eh, Eelec_f_Eh):
    Ei_rel = (level_i.energy_cm + ZPE_i_cm) / EH_TO_CM
    Ef_rel = (level_f.energy_cm + ZPE_f_cm) / EH_TO_CM

    Ei_abs = Eelec_i_Eh + Ei_rel
    Ef_abs = Eelec_f_Eh + Ef_rel

    return Ei_rel, Ef_rel, Ef_abs - Ei_abs


def precompute_dyson_values(dyson_splines, r_use):
    dyson_values = {}

    for key, spline in dyson_splines.items():
        vals = spline(r_use)

        if np.any(~np.isfinite(vals)):
            raise ValueError(
                f"Dyson spline returned non-finite values for key {key}. "
                f"R range used: {r_use.min()} - {r_use.max()}"
            )

        dyson_values[key] = vals

    return dyson_values


def transition_amplitude_component_contraction(
    level_i,
    level_f,
    chi_i,
    chi_f,
    components_i,
    components_f,
    dyson_values_by_key,
    K_values,
    missing_dyson="skip",
):
    matrix_element = 0.0 + 0.0j
    bk_no_rotation = 0.0 + 0.0j
    rot_abs_sum = 0.0
    n_terms = 0
    missing_keys = set()

    for comp_i in components_i:
        for comp_f in components_f:
            keydyson = (comp_i.Omega, comp_f.Omega)

            if keydyson not in dyson_values_by_key:
                missing_keys.add(keydyson)
                continue

            vib_part = bra_ket_weighted(
                vib_ini=chi_i,
                vib_fin=chi_f,
                dyson_values=dyson_values_by_key[keydyson],
            )

            coeff_factor = np.conj(comp_f.coeff) * comp_i.coeff

            rot = rotational_factor(
                J_i=level_i.J,
                J_f=level_f.J,
                Omega_i=comp_i.Omega,
                Omega_f=comp_f.Omega,
                K_values=K_values,
            )

            if np.isclose(rot, 0.0):
                continue

            bk_component = coeff_factor * vib_part
            matrix_element += rot * bk_component
            bk_no_rotation += bk_component
            rot_abs_sum += abs(rot)
            n_terms += 1

    if missing_keys and missing_dyson == "error":
        raise KeyError(
            f"Missing Dyson keys for transition "
            f"i={level_i.duo_index}, f={level_f.duo_index}: "
            f"{sorted(missing_keys)}"
        )

    return matrix_element, bk_no_rotation, rot_abs_sum, n_terms


def build_transition_table(
    levels_neutral,
    levels_cation,
    vib_neutral,
    vib_cation,
    coeff_neutral,
    coeff_cation,
    dyson_splines,
    rvals,
    mask,
    ZPE_neutral_cm,
    ZPE_cation_cm,
    Eelec_neutral_Eh,
    Eelec_cation_Eh,
    K_values,
    missing_dyson="skip",
    missing_coeff="error",
    min_intensity=0.0,
):
    rows = []

    r_use = rvals[mask]
    vib_neutral_use = vib_neutral[mask, :]
    vib_cation_use = vib_cation[mask, :]
    dyson_values_by_key = precompute_dyson_values(dyson_splines, r_use)

    for level_i in tqdm(levels_neutral, desc="Neutral levels"):
        vi = level_i.v
        chi_i = vib_neutral_use[:, vi]
        key_i = level_state_key(level_i)

        if key_i not in coeff_neutral:
            if missing_coeff == "error":
                raise KeyError(f"Missing neutral coefficient key: {key_i}")
            if missing_coeff == "skip":
                continue
            raise ValueError("missing_coeff must be 'error' or 'skip'")

        components_i = coeff_neutral[key_i]

        for level_f in levels_cation:
            vf = level_f.v
            chi_f = vib_cation_use[:, vf]
            key_f = level_state_key(level_f)

            if key_f not in coeff_cation:
                if missing_coeff == "error":
                    raise KeyError(f"Missing cation coefficient key: {key_f}")
                if missing_coeff == "skip":
                    continue
                raise ValueError("missing_coeff must be 'error' or 'skip'")

            components_f = coeff_cation[key_f]

            Ei_rel_Eh, Ef_rel_Eh, DeltaE_Eh = transition_energy_Eh(
                level_i=level_i,
                level_f=level_f,
                ZPE_i_cm=ZPE_neutral_cm,
                ZPE_f_cm=ZPE_cation_cm,
                Eelec_i_Eh=Eelec_neutral_Eh,
                Eelec_f_Eh=Eelec_cation_Eh,
            )

            if DeltaE_Eh <= 0.0:
                continue

            matrix_element, bk_no_rotation, rot_abs_sum, n_terms = transition_amplitude_component_contraction(
                level_i=level_i,
                level_f=level_f,
                chi_i=chi_i,
                chi_f=chi_f,
                components_i=components_i,
                components_f=components_f,
                dyson_values_by_key=dyson_values_by_key,
                K_values=K_values,
                missing_dyson=missing_dyson,
            )

            if n_terms == 0:
                continue

            intensity = np.abs(matrix_element) ** 2

            if intensity <= min_intensity:
                continue

            rows.append(
                {
                    "v_i": level_i.v,
                    "J_i": level_i.J,
                    "Omega_i": level_i.Omega,
                    "Sigma_i": level_i.Sigma,
                    "Lambda_i": level_i.Lambda,
                    "parity_i": level_i.parity,
                    "index_i": level_i.duo_index,
                    "v_f": level_f.v,
                    "J_f": level_f.J,
                    "Omega_f": level_f.Omega,
                    "Sigma_f": level_f.Sigma,
                    "Lambda_f": level_f.Lambda,
                    "parity_f": level_f.parity,
                    "index_f": level_f.duo_index,
                    "Ei_eV": Ei_rel_Eh * EH_TO_EV,
                    "Ef_eV": Ef_rel_Eh * EH_TO_EV,
                    "DeltaE_eV": DeltaE_Eh * EH_TO_EV,
                    "I_raw": intensity,
                    "matrix_real": np.real(matrix_element),
                    "matrix_imag": np.imag(matrix_element),
                    "matrix_abs": np.abs(matrix_element),
                    "bk_real": np.real(bk_no_rotation),
                    "bk_imag": np.imag(bk_no_rotation),
                    "bk_abs": np.abs(bk_no_rotation),
                    "rot_abs_sum": rot_abs_sum,
                    "n_component_terms": n_terms,
                    "Delta_v": level_f.v - level_i.v,
                    "Delta_J": level_f.J - level_i.J,
                    "Delta_Omega": level_f.Omega - level_i.Omega,
                }
            )

    df = pd.DataFrame(rows)

    if len(df) == 0:
        raise RuntimeError("No transitions were generated.")

    return df


# =============================================================================
# Dump
# =============================================================================

def dump_transition_table(df, filename):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    preferred_cols = [
        "v_i", "J_i", "Omega_i", "Sigma_i", "Lambda_i", "parity_i", "index_i",
        "v_f", "J_f", "Omega_f", "Sigma_f", "Lambda_f", "parity_f", "index_f",
        "Ei_eV", "Ef_eV", "DeltaE_eV",
        "I_raw",
        "matrix_real", "matrix_imag", "matrix_abs",
        "bk_real", "bk_imag", "bk_abs",
        "rot_abs_sum", "n_component_terms",
        "Delta_v", "Delta_J", "Delta_Omega",
    ]

    cols = [col for col in preferred_cols if col in df.columns]
    extra = [col for col in df.columns if col not in cols]
    df_out = df[cols + extra]

    df_out.to_csv(
        filename,
        sep=" ",
        index=False,
        float_format="%.10e",
    )

    print(f"Wrote: {filename}")
    print(f"Rows : {len(df_out)}")


# =============================================================================
# Argument parser and main
# =============================================================================

def build_K_values(DJ_values):
    """
    Expands maximum DJ/K values into the allowed angular ranks.

    Example:
        -DJ 1.5      -> [0.5, 1.5]
        -DJ 2.5      -> [0.5, 1.5, 2.5]
        -DJ 0.5 1.5  -> [0.5, 1.5]
    """

    K_set = set()

    for DJ in DJ_values:
        if DJ < 0.5:
            raise ValueError(f"DJ must be >= 0.5, got {DJ}")

        nsteps = int(round(DJ - 0.5))

        if not np.isclose(0.5 + nsteps, DJ):
            raise ValueError(
                f"DJ={DJ} is not compatible with half-integer sequence "
                "0.5, 1.5, 2.5, ..."
            )

        for n in range(nsteps + 1):
            K_set.add(0.5 + n)

    return sorted(K_set)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build raw rovibronic photoionization transition table."
    )

    parser.add_argument("-ZPE", type=float, nargs=2, required=True,
                        help="Neutral and cation ZPE values in cm^-1.")
    parser.add_argument("-Etot", type=float, nargs=2, required=True,
                        help="Neutral and cation electronic/reference energies in Eh.")
    parser.add_argument("-nVJtot", type=int, nargs=2, required=True,
                        help="Total number of vibrational functions and total J count used for reading.")
    parser.add_argument("-DJ", type=float, nargs="+", required=True,
                        help="Allowed angular K/DJ values, e.g. -DJ 0.5 1.5")
    parser.add_argument("-Nvib", "--numvib", type=int, nargs=2, required=True,
                        help="Number of neutral and cation vibrational states to include.")
    parser.add_argument("-NJ", "--numJ", type=int, nargs=2, default=None,
                        help="Optional J cutoffs for neutral and cation levels.")
    parser.add_argument("-maskR", type=float, nargs=2, required=True,
                        help="R range used in the vibrational overlap.")
    parser.add_argument("-npts", type=int, required=True,
                        help="Number of radial grid points in the vibrational functions.")
    parser.add_argument("-p", "--path", type=str, required=True,
                        help="Base path.")
    parser.add_argument("-dsys", type=str, nargs=2, required=True,
                        help="Neutral and cation system directory/prefix names.")
    parser.add_argument("-dypath", "--pathtody", type=str, required=True,
                        help="Path to Dyson files dyson_XX_YY.dat.")
    parser.add_argument("-stts", "--states", nargs=2, type=str, required=True,
                        help="Neutral and cation state labels, e.g. PHX3Sm PHMX2P.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output transition table.")

    parser.add_argument("--dyson-mode", type=str, default="quadrature",
                        choices=["quadrature", "sum"],
                        help="How to combine Dyson files sharing the same Omega key.")
    parser.add_argument("--spline", type=str, default="cubic",
                        choices=["cubic", "linear"],
                        help="Dyson interpolation type.")
    parser.add_argument("--extrapolate", action="store_true",
                        help="Allow Dyson spline extrapolation outside its R grid.")
    parser.add_argument("--missing-dyson", type=str, default="skip",
                        choices=["skip", "error"],
                        help="What to do if a Dyson key is missing.")
    parser.add_argument("--missing-coeff", type=str, default="error",
                        choices=["skip", "error"],
                        help="What to do if a rovibronic coefficient key is missing.")
    parser.add_argument("--min-intensity", type=float, default=0.0,
                        help="Minimum raw intensity to write.")

    return parser.parse_args()


def main():
    args = parse_args()

    base_path = Path(args.path)
    dyson_path = Path(args.pathtody)

    system_neutral, system_cation = args.dsys
    state_neutral, state_cation = args.states

    nvib_total, _ = args.nVJtot
    nvib_neutral, nvib_cation = args.numvib

    jmax_neutral = None
    jmax_cation = None
    if args.numJ is not None:
        jmax_neutral, jmax_cation = args.numJ

    ZPE_neutral_cm, ZPE_cation_cm = args.ZPE
    Eelec_neutral_Eh, Eelec_cation_Eh = args.Etot
    rmin, rmax = args.maskR
    K_values = build_K_values(args.DJ)

    dyson_index_to_omega = build_dyson_index_to_omega()

    if state_neutral not in dyson_index_to_omega:
        raise KeyError(
            f"Unknown neutral state {state_neutral}. Available: {list(dyson_index_to_omega)}"
        )

    if state_cation not in dyson_index_to_omega:
        raise KeyError(
            f"Unknown cation state {state_cation}. Available: {list(dyson_index_to_omega)}"
        )

    if args.output is None:
        k_tag = "_K" + "-".join(f"{k:g}" for k in K_values)
        output_file = base_path / f"{state_neutral}_{state_cation}{k_tag}.dat"
    else:
        output_file = Path(args.output)

    print("----------------------------------------------------------------")
    print("Reading vibrational eigenfunctions")

    vib_neutral_file = resolve_system_file(base_path, system_neutral, "vibeigenvect_vib.chk")
    vib_cation_file = resolve_system_file(base_path, system_cation, "vibeigenvect_vib.chk")

    vib_neutral = parse_duo_vib_einfun(vib_neutral_file, args.npts, nvib_total)
    vib_cation = parse_duo_vib_einfun(vib_cation_file, args.npts, nvib_total)

    vib_neutral = vib_neutral[:, :nvib_neutral]
    vib_cation = vib_cation[:, :nvib_cation]

    print(f"Neutral vib shape: {vib_neutral.shape}")
    print(f"Cation  vib shape: {vib_cation.shape}")

    print("----------------------------------------------------------------")
    print("Reading DUO rovibronic levels")

    levels_neutral_file = resolve_system_file(base_path, system_neutral, "rovibronic_energies.dat")
    levels_cation_file = resolve_system_file(base_path, system_cation, "rovibronic_energies.dat")

    levels_neutral = read_duo_levels(
        levels_neutral_file,
        system=system_neutral,
        nvib_max=nvib_neutral,
        j_max=jmax_neutral,
    )

    levels_cation = read_duo_levels(
        levels_cation_file,
        system=system_cation,
        nvib_max=nvib_cation,
        j_max=jmax_cation,
    )

    print(f"Neutral levels: {len(levels_neutral)}")
    print(f"Cation  levels: {len(levels_cation)}")

    if len(levels_neutral) == 0:
        raise RuntimeError("No neutral DUO levels read.")

    if len(levels_cation) == 0:
        raise RuntimeError("No cation DUO levels read.")

    print("----------------------------------------------------------------")
    print("Reading DUO coefficients")

    coeff_neutral_file = resolve_system_file(base_path, system_neutral, "vibeigenvect_vectors.chk")
    coeff_cation_file = resolve_system_file(base_path, system_cation, "vibeigenvect_vectors.chk")

    coeff_neutral = read_coefficients(coeff_neutral_file)
    coeff_cation = read_coefficients(coeff_cation_file)

    print(f"Neutral coefficient keys: {len(coeff_neutral)}")
    print(f"Cation  coefficient keys: {len(coeff_cation)}")

    print("----------------------------------------------------------------")
    print("Reading R grid")

    rgrid_file = base_path / "Dipole_moment_functions.dat"
    if not rgrid_file.exists():
        raise FileNotFoundError(f"Missing R-grid file: {rgrid_file}")

    rvals = np.loadtxt(rgrid_file, usecols=0)
    mask = (rvals > rmin) & (rvals < rmax)

    if np.sum(mask) == 0:
        raise RuntimeError(
            f"Empty R mask. Requested range: {rmin} - {rmax}. "
            f"Available range: {np.min(rvals)} - {np.max(rvals)}"
        )

    print(f"R used: {rvals[mask].min()} - {rvals[mask].max()}")
    print(f"R points used: {np.sum(mask)}")

    print("----------------------------------------------------------------")
    print("Reading Dyson splines")

    dyson_splines, _, _ = read_dyson_splines(
        d_neutral=dyson_index_to_omega[state_neutral],
        d_cation=dyson_index_to_omega[state_cation],
        dyson_path=dyson_path,
        combine_mode=args.dyson_mode,
        spline_type=args.spline,
        bohr_to_angstrom=True,
        extrapolate=args.extrapolate,
    )

    print(f"Dyson spline keys: {len(dyson_splines)}")

    print("----------------------------------------------------------------")
    print("Building transition table")

    df = build_transition_table(
        levels_neutral=levels_neutral,
        levels_cation=levels_cation,
        vib_neutral=vib_neutral,
        vib_cation=vib_cation,
        coeff_neutral=coeff_neutral,
        coeff_cation=coeff_cation,
        dyson_splines=dyson_splines,
        rvals=rvals,
        mask=mask,
        ZPE_neutral_cm=ZPE_neutral_cm,
        ZPE_cation_cm=ZPE_cation_cm,
        Eelec_neutral_Eh=Eelec_neutral_Eh,
        Eelec_cation_Eh=Eelec_cation_Eh,
        K_values=K_values,
        missing_dyson=args.missing_dyson,
        missing_coeff=args.missing_coeff,
        min_intensity=args.min_intensity,
    )

    if np.any(df["index_i"] <= 0):
        raise RuntimeError("Found neutral DUO index <= 0 in output table.")

    if np.any(df["index_f"] <= 0):
        raise RuntimeError("Found cation DUO index <= 0 in output table.")

    print(f"Transitions: {len(df)}")
    print(f"Total I_raw: {df['I_raw'].sum():.10e}")
    print(f"Max   I_raw: {df['I_raw'].max():.10e}")

    print("----------------------------------------------------------------")
    print("Writing output")

    dump_transition_table(df, output_file)

    print("----------------------------------------------------------------")
    print("Done")


if __name__ == "__main__":
    main()
