# Coolant inputs and corrected Gnielinski heat transfer. All properties use SI.

import math


# Select one coolant or several to compare in the same run.
SELECTED_COOLANTS = ["IPA", "Ethanol"]  # Example: ["IPA", "Ethanol"]
MASS_FLOW = 1.456 * 0.295  # kg/s, shared by the comparison cases
INLET_TEMPERATURE = 293.15  # K, entering at the nozzle exit


def ipa_density(temperature):
    # Kerscher et al., doi:10.1007/s10765-023-03294-z, Eq. 3, Table 2.
    return (1161.73 - 2.5804 * temperature + 7.2987e-3 * temperature**2
            - 9.8196e-6 * temperature**3)


def ipa_viscosity(temperature):
    # Same source, Eq. 5, Table 6. Convert mPa s to Pa s.
    log_viscosity = (math.log(1.23842e-16) + 4.37443e4 / temperature
                     - 2.17822e7 / temperature**2 + 5.03956e9 / temperature**3
                     - 4.30435e11 / temperature**4)
    return math.exp(log_viscosity) / 1000


# Installed RPA 2.3 properties.inp, C2H5OH(L), interpolated to 3.447 MPa.
# Columns: temperature (K), density (kg/m^3), viscosity (Pa s).
ETHANOL_TABLE = [
    (253, 825.033536, 0.00296129634),
    (273, 808.854614, 0.00185462742),
    (293, 792.655056, 0.00122371628),
    (313, 775.899544, 0.00084310873),
    (333, 758.147654, 0.000602201076),
    (353, 739.004068, 0.00044312966),
    (373, 718.082620, 0.000333941266),
    (393, 694.942038, 0.000256243968),
    (413, 668.939778, 0.000199050098),
]


def ethanol_property(temperature, column):
    for lower, upper in zip(ETHANOL_TABLE, ETHANOL_TABLE[1:]):
        if lower[0] <= temperature <= upper[0]:
            fraction = (temperature - lower[0]) / (upper[0] - lower[0])
            return lower[column] + fraction * (upper[column] - lower[column])
    raise ValueError("Ethanol table requires 253-413 K")


# Copy a block to add a coolant. Density and viscosity accept a number or a
# function of temperature. Heat capacity and conductivity are constant inputs.
# Use measured properties and their valid range; do not reuse IPA fits for another fuel.
COOLANTS = {
    "IPA": {
        "density": ipa_density,  # kg/m^3
        "viscosity": ipa_viscosity,  # Pa s
        "specific_heat": 2600.0,  # J/(kg K)
        "conductivity": 0.14,  # W/(m K)
        "temperature_range": (278, 403),  # K, bulk-property range
        "wall_property_limit": 403.0,  # K; None rejects walls outside the range
    },
    "Ethanol": {
        "density": lambda temperature: ethanol_property(temperature, 1),
        "viscosity": lambda temperature: ethanol_property(temperature, 2),
        # Inlet values from the same RPA table, held constant in the heat balance.
        "specific_heat": 2506.9879039,
        "conductivity": 0.16816263528,
        "temperature_range": (253, 413),
        "wall_property_limit": 413.0,
    },
}


def coolant_inputs(name):
    fuel = COOLANTS[name]
    return {
        "coolant_name": name,
        "coolant_mass_flow": MASS_FLOW,
        "coolant_temperature": INLET_TEMPERATURE,
        **{f"coolant_{key}": value for key, value in fuel.items()},
    }


def liquid_properties(temperature, inputs):
    lower, upper = inputs["coolant_temperature_range"]
    if not lower <= temperature <= upper:
        raise ValueError(f"{inputs['coolant_name']} properties require {lower}-{upper} K")
    properties = {}
    for name in ("density", "viscosity", "specific_heat", "conductivity"):
        value = inputs[f"coolant_{name}"]
        value = value(temperature) if callable(value) else value
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"Coolant {name} must be finite and positive")
        properties[f"coolant_{name}"] = value
    return properties


def channel_flow(channel, inputs, temperature, wall_temperature):
    properties = liquid_properties(temperature, inputs)
    mass_velocity = inputs["coolant_mass_flow"] / (
        inputs["number_of_channels"] * channel["flow_area"]
    )
    reynolds = mass_velocity * channel["hydraulic_diameter"] / properties["coolant_viscosity"]
    prandtl = (properties["coolant_specific_heat"] * properties["coolant_viscosity"]
               / properties["coolant_conductivity"])
    if not 2300 <= reynolds <= 5e6 or not 0.5 <= prandtl <= 2000:
        raise ValueError("Gnielinski requires 2300 <= Re <= 5e6 and 0.5 <= Pr <= 2000")

    # Limit property evaluation only; the solved wall temperature is never clipped.
    limit = inputs["coolant_wall_property_limit"]
    wall_property_temperature = wall_temperature
    if limit is not None:
        liquid_properties(limit, inputs)
        wall_property_temperature = min(wall_temperature, limit)
    wall_properties = liquid_properties(wall_property_temperature, inputs)
    wall_prandtl = (wall_properties["coolant_specific_heat"] * wall_properties["coolant_viscosity"]
                   / wall_properties["coolant_conductivity"])
    correction = (prandtl / wall_prandtl)**0.11

    friction = (0.79 * math.log(reynolds) - 1.64)**-2  # Smooth-channel Darcy factor
    nusselt = (friction / 8) * (reynolds - 1000) * prandtl / (
        1 + 12.7 * math.sqrt(friction / 8) * (prandtl**(2 / 3) - 1)
    ) * correction
    coefficient = nusselt * properties["coolant_conductivity"] / channel["hydraulic_diameter"]
    return {
        **properties,
        "coolant_name": inputs["coolant_name"],
        "coolant_velocity": mass_velocity / properties["coolant_density"],
        "reynolds_number": reynolds,
        "prandtl_number": prandtl,
        "nusselt_number": nusselt,
        "coolant_coefficient": coefficient,
        "hot_surface_correction": correction,
        "wall_property_temperature": wall_property_temperature,
        "wall_properties_limited": wall_property_temperature < wall_temperature,
    }


def sidewall_cooling(channel, coefficient, wall_conductivity):
    # Straight rectangular ribs with insulated tips; include floor and sidewalls.
    fin_parameter = channel["height"] * math.sqrt(
        2 * coefficient / (wall_conductivity * channel["rib_width"])
    )
    efficiency = math.tanh(fin_parameter) / fin_parameter
    area_ratio = (channel["width"] + 2 * channel["height"] * efficiency) / channel["hot_wall_pitch"]
    return {
        "fin_efficiency": efficiency,
        "effective_cooled_area_ratio": area_ratio,
        "effective_coolant_coefficient": coefficient * area_ratio,
    }
