# Engine and channel inputs. Coolant selection and properties are in Coolant.py.
# Channel dimensions and wall thickness are in mm; other inputs use SI units.

import csv
from pathlib import Path

import Contour
import Coolant
import Wall


def model_inputs():
    return {
        "gas_mass_flow": 1.456,  # kg/s
        # Chamber thermodynamics from the archived IPA/N2O CEA report.
        "stagnation_temperature": 2316.97,
        "gas_specific_heat": 1953.4,
        "gamma": 1.2592,
        # Transport assumptions retained from v1; absent from the CEA report.
        "gas_viscosity": 14.7e-6,
        "gas_conductivity": 0.01662,
        "gas_viscosity_exponent": 0.6,
        "gas_recovery_factor": 0.90,  # Engineering approximation
        "number_of_channels": 40,
        "channel_height": 1.5,
        "chamber_channel_width": 4.0,
        "throat_channel_width": 1.5,
        "exit_channel_width": 2.0,
        "wall_thickness": 1.0,
        "wall_conductivity": 344.0,
    }


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    if not Coolant.SELECTED_COOLANTS:
        raise ValueError("Select at least one coolant in Coolant.py")
    contour = Contour.load_contour()
    output = Path(__file__).parent / "outputs"
    output.mkdir(exist_ok=True)
    summaries = []

    for name in Coolant.SELECTED_COOLANTS:
        inputs = {**model_inputs(), **Coolant.coolant_inputs(name)}
        profile = Wall.wall_temperature_profile(inputs, contour, station_count=400)
        save_csv(output / f"profile_{name}.csv", profile)
        peak = max(profile, key=lambda point: point["hot_wall_temperature"])
        summary = {
            "coolant": name,
            "peak_wall_K": peak["hot_wall_temperature"],
            "peak_x_mm": peak["x"],
            "coolant_outlet_K": profile[0]["coolant_outlet_temperature"],
            "heat_kW": sum(point["heat_rate"] for point in profile) / 1000,
            "wall_property_limited_segments": sum(point["wall_properties_limited"] for point in profile),
            "segments_below_Re_10000": sum(point["reynolds_number"] < 10000 for point in profile),
        }
        summaries.append(summary)
        print(f"{name}: peak wall {summary['peak_wall_K']:.2f} K at {summary['peak_x_mm']:.2f} mm, "
              f"coolant outlet {summary['coolant_outlet_K']:.2f} K, heat {summary['heat_kW']:.2f} kW")
        print(f"  Wall-property limit used in {summary['wall_property_limited_segments']}/{len(profile)} segments; "
              f"Re < 10000 in {summary['segments_below_Re_10000']}/{len(profile)}.")

    if summaries:
        save_csv(output / "summary.csv", summaries)
    print(f"Saved results in {output}")


if __name__ == "__main__":
    main()
