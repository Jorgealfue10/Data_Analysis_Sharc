import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Argumentos de línea de comandos
parser = argparse.ArgumentParser(description="Plot photoionization spectrum with peak labels avoiding overlap")
parser.add_argument("filename", help="Archivo de entrada con el espectro")
parser.add_argument("--output", help="Archivo de salida opcional para guardar la figura", default=None)
parser.add_argument("--min_distance", type=float, default=0.1, help="Distancia mínima entre etiquetas de picos (eV)")

args = parser.parse_args()

# Leer el archivo y limpiar las líneas
with open("spectrum.out", "r") as f:
    lines = f.readlines()

# Separar los valores asegurando consistencia de columnas
data = [line.strip().split() for line in lines]
lengths = [len(row) for row in data]
most_common_length = max(set(lengths), key=lengths.count)
filtered = [row for row in data if len(row) == most_common_length]

# Convertir a DataFrame con valores float
df = pd.DataFrame(filtered).apply(pd.to_numeric, errors="coerce").dropna()

# Separar columnas
energy = df.iloc[:, 0]
states = df.iloc[:, 1:-1]
total = df.iloc[:, -1]

# Graficar
plt.figure(figsize=(16, 12))
# for i in range(states.shape[1]):
#     plt.fill_between(energy, states.iloc[:, i], alpha=0.4, label=f"State {i+1}")

# for i in range(states.shape[1]):
#     y = states.iloc[:, i]
#     plt.fill_between(energy, y, alpha=0.4, label=f"State {i+1}")

#     # Detectar pico y etiquetar
#     idx_max = y.idxmax()
#     if y[idx_max] > 0.01 * total.max():  # Solo etiquetar picos significativos (>1% del máximo total)
#         x_peak = energy[idx_max]
#         y_peak = y[idx_max]
#         plt.text(x_peak, y_peak, f"{x_peak:.2f}", ha="center", va="bottom", fontsize=8, rotation=45)

labeled_peaks = []  # lista de posiciones etiquetadas

for i in range(states.shape[1]):
    y = states.iloc[:, i]
    plt.fill_between(energy, y, alpha=0.4, label=f"State {i+1}")

    idx_max = y.idxmax()
    x_peak = energy[idx_max]
    y_peak = y[idx_max]

    # Condición: pico significativo y no demasiado cerca de etiquetas previas
    if y_peak > 0.01 * total.max() and all(abs(x_peak - xp) > args.min_distance for xp in labeled_peaks):
        labeled_peaks.append(x_peak)
        plt.text(x_peak,total.max(), f"{x_peak:.2f}", ha="center", va="bottom", fontsize=14, rotation=45)
        plt.vlines(x_peak, 0, total.max(), color="red", linestyle="--", linewidth=0.5)  # Línea vertical para el pico


plt.plot(energy, total, color="black", linewidth=2, label="Total")

plt.xlabel("Energy (eV)")
plt.ylabel("Intensity (a.u.)")
plt.title(r"PH$_3\rightarrow$PH$_3^+$ Photoionization Spectrum")
# plt.legend(loc="upper right", fontsize="small", ncol=2)
plt.xlim(8,20)
plt.grid(True)
plt.tight_layout()

if args.output:
    plt.savefig(args.output)
    plt.show()
else:
    plt.show()
