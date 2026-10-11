# Plot the saved IPA, ethanol, and RPA hot-wall temperature distributions.
# Run main.py with both coolants selected, then run this file.

import csv
from pathlib import Path

import matplotlib.pyplot as plt

import Contour


def main():
    root = Path(__file__).parent
    output = root / "outputs"
    curves = []
    for name, color in (("IPA", "#0072b2"), ("Ethanol", "#d55e00")):
        with (output / f"profile_{name}.csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        curves.append((
            name, color,
            [float(row["x"]) for row in rows],
            [float(row["hot_wall_temperature"]) for row in rows],
        ))

    rpa = []
    for line in (root / "v1" / "RPA output.txt").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 13:
            continue
        try:
            values = [float(value) for value in fields[:13]]
        except ValueError:
            continue
        rpa.append((values[0], values[6]))  # Axial position (mm), gas-side wall (K)
    if not rpa:
        raise ValueError("No thermal stations found in the saved RPA output")

    contour = Contour.load_contour()
    throat_x = contour.throat()["x"]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    for axis in axes:
        for name, color, x, temperature in curves:
            axis.plot(x, temperature, color=color, linewidth=2, label=f"{name} model")
        axis.plot([row[0] for row in rpa], [row[1] for row in rpa],
                  "k.--", linewidth=1.3, markersize=4, label="RPA — ethanol proxy")
        axis.axvline(throat_x, color="#888888", linestyle=":", linewidth=1)
        axis.set_xlabel("Axial position from injector (mm)")
        axis.set_ylabel("Gas-side wall temperature (K)")
        axis.grid(alpha=0.2)
    axes[0].set(title="Full chamber and nozzle", xlim=(contour.x_min, contour.x_max))
    axes[1].set(title="Throat to nozzle exit", xlim=(throat_x, contour.x_max))
    nozzle_temperatures = [temperature for _, _, x, temperatures in curves
                           for position, temperature in zip(x, temperatures) if position >= throat_x]
    nozzle_temperatures.extend(temperature for x, temperature in rpa if x >= throat_x)
    axes[1].set_ylim(min(nozzle_temperatures) - 20, max(nozzle_temperatures) + 20)
    axes[0].legend(fontsize=9)
    figure.suptitle("Osiris hot-wall temperature distributions", fontsize=15)
    figure.text(0.5, 0.02,
                "Dotted vertical line: throat. Models use corrected Gnielinski; RPA includes film cooling and radiation.",
                ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.07, 1, 0.94))
    path = output / "temperature_distributions.png"
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
