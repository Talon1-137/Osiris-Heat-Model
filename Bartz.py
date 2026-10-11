# Gas-side heat transfer from Rocket Propulsion Elements, Eq. 8-22.
# All inputs use SI units. Gas properties are referenced to the chamber.


def heat_transfer_coefficient(wall_temperature, gas, inputs):
    mean_temperature = (gas["static_temperature"] + wall_temperature) / 2

    # Ideal-gas density at the same pressure; viscosity follows mu ~ T**omega.
    density_ratio = gas["static_temperature"] / mean_temperature
    viscosity_ratio = (
        mean_temperature / inputs["stagnation_temperature"]
    )**inputs["gas_viscosity_exponent"]
    prandtl = (
        inputs["gas_specific_heat"] * inputs["gas_viscosity"]
        / inputs["gas_conductivity"]
    )

    return (
        0.026 / gas["diameter"]**0.2
        * inputs["gas_specific_heat"] * inputs["gas_viscosity"]**0.2
        / prandtl**0.6
        * gas["mass_velocity"]**0.8
        * density_ratio * viscosity_ratio**0.2
    )
