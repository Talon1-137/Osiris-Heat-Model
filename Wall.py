# Steady radial heat flow, Rocket Propulsion Elements, Eqs. 8-15 to 8-19.
# Each segment heats the coolant according to Eq. 8-25.

import math

import Bartz
import Coolant
import Geometry
import Thermo


def heat_balance(gas_temperature, coolant_inlet_temperature, gas_coefficient,
                 coolant_coefficient, wall_thickness, wall_conductivity,
                 coolant_warming_resistance=0):
    # Thickness is in m. All resistances are per unit hot-wall area, m^2 K/W.
    wall_resistance = wall_thickness / wall_conductivity
    heat_flux = (gas_temperature - coolant_inlet_temperature) / (
        1 / gas_coefficient + wall_resistance + 1 / coolant_coefficient
        + coolant_warming_resistance
    )
    coolant_bulk = (
        coolant_inlet_temperature + heat_flux * coolant_warming_resistance
    )
    return {
        "heat_flux": heat_flux,
        "hot_wall_temperature": gas_temperature - heat_flux / gas_coefficient,
        "coolant_wall_temperature": coolant_bulk + heat_flux / coolant_coefficient,
        "coolant_bulk_temperature": coolant_bulk,
    }


def solve_wall(inputs, gas, coolant_inlet, effective_coefficient, warming_resistance):
    lower = min(coolant_inlet, gas["recovery_temperature"])
    upper = max(coolant_inlet, gas["recovery_temperature"])

    def evaluate(wall_temperature):
        coefficient = Bartz.heat_transfer_coefficient(wall_temperature, gas, inputs)
        result = heat_balance(
            gas["recovery_temperature"], coolant_inlet, coefficient,
            effective_coefficient, inputs["wall_thickness"] / 1000,
            inputs["wall_conductivity"], warming_resistance,
        )
        result["gas_coefficient"] = coefficient
        return result

    for _ in range(80):
        temperature = (lower + upper) / 2
        result = evaluate(temperature)
        if result["hot_wall_temperature"] > temperature:
            lower = temperature
        else:
            upper = temperature
        if upper - lower < 1e-7:
            break
    else:
        raise RuntimeError("Gas-side wall-temperature iteration did not converge")

    return evaluate((lower + upper) / 2)


def point_wall_temperature(inputs, gas, channel, heated_area, coolant_inlet):
    heat_capacity_rate = inputs["coolant_mass_flow"] * inputs["coolant_specific_heat"]
    warming_resistance = heated_area / (2 * heat_capacity_rate)
    bulk_temperature = coolant_inlet
    wall_temperature = coolant_inlet

    # Couple bulk properties and the coolant-facing wall correction to the balance.
    for _ in range(80):
        flow = Coolant.channel_flow(channel, inputs, bulk_temperature, wall_temperature)
        cooling = Coolant.sidewall_cooling(
            channel, flow["coolant_coefficient"], inputs["wall_conductivity"],
        )
        result = solve_wall(
            inputs, gas, coolant_inlet, cooling["effective_coolant_coefficient"],
            warming_resistance,
        )
        if (
            abs(result["coolant_bulk_temperature"] - bulk_temperature) < 1e-7
            and abs(result["coolant_wall_temperature"] - wall_temperature) < 1e-7
        ):
            break
        bulk_temperature = result["coolant_bulk_temperature"]
        wall_temperature = result["coolant_wall_temperature"]
    else:
        raise RuntimeError("Coolant bulk/wall-temperature iteration did not converge")

    heat_rate = result["heat_flux"] * heated_area
    coolant_outlet = coolant_inlet + heat_rate / heat_capacity_rate
    # Reject extrapolation even if only the very end of a segment exceeds the range.
    Coolant.liquid_properties(coolant_outlet, inputs)
    result.update(flow)
    result.update(cooling)
    result.update({
        "heated_area": heated_area,
        "heat_rate": heat_rate,
        "coolant_inlet_temperature": coolant_inlet,
        "coolant_outlet_temperature": coolant_outlet,
    })
    return result


def wall_temperature_profile(inputs, contour, station_count=100):
    positive_inputs = (
        "gas_mass_flow", "stagnation_temperature", "gas_specific_heat",
        "gas_viscosity", "gas_conductivity", "gamma", "number_of_channels",
        "channel_height", "chamber_channel_width", "throat_channel_width",
        "exit_channel_width", "wall_thickness", "wall_conductivity",
        "coolant_mass_flow", "coolant_temperature",
        "coolant_conductivity", "coolant_specific_heat",
    )
    for name in positive_inputs:
        if not math.isfinite(inputs[name]) or inputs[name] <= 0:
            raise ValueError(f"{name} must be finite and positive")
    channels = inputs["number_of_channels"]
    if isinstance(channels, bool) or int(channels) != channels or inputs["gamma"] <= 1:
        raise ValueError("Use a whole channel count and gamma greater than one")
    if not math.isfinite(inputs["gas_viscosity_exponent"]):
        raise ValueError("gas_viscosity_exponent must be finite")

    profile = []
    coolant_temperature = inputs["coolant_temperature"]
    for station in reversed(Geometry.axial_stations(contour, station_count)):
        x = station["x"]
        gas = Thermo.local_gas(x, contour, inputs)
        channel = Geometry.channel_geometry(x, contour, inputs)
        result = point_wall_temperature(
            inputs, gas, channel, station["heated_area"], coolant_temperature
        )
        result.update({
            "x": x,
            "radius_mm": contour.radius_at(x),
            "area_ratio": gas["area_ratio"],
            "mach": gas["mach"],
            "gas_static_temperature": gas["static_temperature"],
            "adiabatic_wall_temperature": gas["recovery_temperature"],
        })
        profile.append(result)
        coolant_temperature = result["coolant_outlet_temperature"]

    return profile[::-1]
