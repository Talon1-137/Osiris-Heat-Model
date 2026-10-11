# Constant-gamma, isentropic gas flow using the area read from the DXF.

import math


def area_ratio(mach, gamma):
    factor = 2 / (gamma + 1) * (1 + (gamma - 1) * mach**2 / 2)
    exponent = (gamma + 1) / (2 * (gamma - 1))
    return factor**exponent / mach


def mach_from_area_ratio(ratio, gamma, supersonic):
    if ratio < 1 or gamma <= 1:
        raise ValueError("Area ratio must be at least one and gamma must exceed one")
    if math.isclose(ratio, 1, rel_tol=0, abs_tol=1e-10):
        return 1.0

    if supersonic:
        low, high = 1.0, 2.0
        while area_ratio(high, gamma) < ratio:
            high *= 2
    else:
        low, high = 1e-10, 1.0

    for _ in range(80):
        mach = (low + high) / 2
        below_target = area_ratio(mach, gamma) < ratio
        if below_target == supersonic:
            low = mach
        else:
            high = mach
    return (low + high) / 2


def local_gas(x, contour, inputs):
    throat = contour.throat()
    area = contour.area_at(x)
    ratio = area / throat["area"]
    gamma = inputs["gamma"]
    mach = mach_from_area_ratio(ratio, gamma, x > throat["x"])
    mach_term = 1 + (gamma - 1) * mach**2 / 2
    static_temperature = inputs["stagnation_temperature"] / mach_term
    # Explicit engineering estimate, independent of unverified transport inputs.
    # For this combustion-gas model, 0 <= r <= 1; no temperature clipping.
    recovery_factor = inputs["gas_recovery_factor"]
    if not math.isfinite(recovery_factor) or not 0 <= recovery_factor <= 1:
        raise ValueError("gas_recovery_factor must be between zero and one for this model")
    recovery_temperature = static_temperature + recovery_factor * (
        inputs["stagnation_temperature"] - static_temperature
    )
    return {
        "area_ratio": ratio,
        "mach": mach,
        "diameter": 2 * contour.radius_at(x) / 1000,
        "mass_velocity": inputs["gas_mass_flow"] / area,
        "static_temperature": static_temperature,
        "recovery_temperature": recovery_temperature,
    }
