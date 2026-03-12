import numpy as np
import matplotlib.pyplot as plt
import os ; import re
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
import argparse

def convolve_spectrum(evals,relInt,dE,sigma_eV,emin=None,emax=None,normalize=False):

    energies = evals.flatten()
    intensities = relInt.flatten()

    if emin is None: emin = np.min(energies)
    if emax is None: emax = np.max(energies)

    E_grid = np.arange(emin, emax, dE)
    spectrum = np.zeros_like(E_grid)

    for E_k, I_k in zip(energies, intensities):
        idx = np.abs(E_grid - E_k) <= dE
        spectrum[idx] += I_k

    sigma_bins = sigma_eV / dE
    bins = np.append(E_grid, E_grid[-1] + dE)

    # spectrum, _ = np.histogram(
    #     energies, 
    #     bins=bins, 
    #     weights=intensities
    # )

    kernel = Gaussian1DKernel(sigma_bins)
    spectrum_conv = convolve(spectrum, kernel, normalize_kernel=True, boundary='extend')

    if normalize:
        spectrum_conv = spectrum_conv / np.max(spectrum_conv)
        spectrum = spectrum / np.max(spectrum)

    return E_grid, spectrum, spectrum_conv

def upper_envelope(E, spectrum):

    peaks, _ = find_peaks(spectrum)

    E_peaks = E[peaks]
    I_peaks = spectrum[peaks]

    f = interp1d(E_peaks, I_peaks, kind='cubic', fill_value="extrapolate")

    envelope = f(E)

    return envelope


def lorentzian(E, E0, gamma, I):
    return I * (gamma**2 / ((E - E0)**2 + gamma**2))

def spectrum_lorentzian(evals, relInt, dE, gamma, emin=None, emax=None, normalize=False):

    energies = evals.flatten()
    intensities = relInt.flatten()

    if emin is None:
        emin = np.min(energies)

    if emax is None:
        emax = np.max(energies)

    # grilla de energía
    E_grid = np.arange(emin, emax, dE)

    # guardar cada lorentziana
    lorentzians = []

    for Ek, Ik in zip(energies, intensities):

        L = lorentzian(E_grid, Ek, gamma, Ik)
        lorentzians.append(L)

    lorentzians = np.array(lorentzians)

    # sumar todas
    spectrum = np.sum(lorentzians, axis=0)

    if normalize:
        spectrum /= np.max(spectrum)

    return E_grid, spectrum, lorentzians

def main():
    parser = argparse.ArgumentParser(
        description='Script to convolve a spectrum'
    )

    parser.add_argument('-i', '--input', type=str, help='Input file')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('--emin', type=float, default=None, help='Minimum energy')
    parser.add_argument('--emax', type=float, default=None, help='Maximum energy')
    parser.add_argument('--dE', type=float, default=0.01, help='Energy step')
    parser.add_argument('--sigma_eV', type=float, default=0.1, help='Sigma in eV')
    parser.add_argument('-N', '--normalize', action='store_true', help='Normalize spectrum')

    args = parser.parse_args()

    fileinp = args.input
    fileout = args.output
    emin = args.emin
    emax = args.emax
    dE = args.dE
    sigma_eV = args.sigma_eV
    if args.normalize:
        normalize = True
    else:
        normalize = False

    data = np.loadtxt(fileinp, skiprows=1)
    evals, relInt = data[:,-2], data[:,-1]
    # E_grid, spectrum, spectrum_conv = convolve_spectrum(evals,relInt,dE,sigma_eV,emin,emax,normalize)
    # E_grid, spectrum, spectrum_conv = convolve_spectrum(evals,relInt,dE,sigma_eV,emin,emax,normalize)
    E, spec, lorentzians = spectrum_lorentzian(evals, relInt, dE, sigma_eV, emin, emax, normalize)
    # np.savetxt(fileout, np.column_stack((E_grid, envelope, spectrum, spectrum_conv)))
    np.savetxt(fileout, np.column_stack((E,spec)))


if __name__ == '__main__':
    main()