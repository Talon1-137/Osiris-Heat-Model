# Osiris heat model

```powershell
python -m pip install -r requirements.txt
python main.py
```

Edit engine and channel inputs in `main.py`. Choose coolants and edit their
properties, mass flow, and inlet temperature in `Coolant.py`.

## Compare coolants

At the top of `Coolant.py`, select one or more names from `COOLANTS`:

```python
SELECTED_COOLANTS = ["IPA", "Ethanol"]
```

Run `main.py`. Each coolant gets `outputs/profile_<name>.csv`, and
`outputs/summary.csv` compares peak wall temperature, coolant outlet
temperature, and absorbed heat. Only the selected cases are run; the summary
is replaced each time. Previously generated individual profiles remain.

To add a coolant, copy a `COOLANTS` block and enter density (kg/mÂ³), viscosity
(Pa s), specific heat (J/(kg K)), conductivity (W/(m K)), and the valid
temperature range (K). Density and viscosity can be numbers or functions of
temperature. Specific heat and conductivity are constant inputs.

```python
"MyFuel": {
    "density": density_in_kg_m3,
    "viscosity": viscosity_in_Pa_s,
    "specific_heat": cp_in_J_kgK,
    "conductivity": conductivity_in_W_mK,
    "temperature_range": (minimum_K, maximum_K),
    "wall_property_limit": None,
},
```

Use actual property values in place of these placeholders. Constant viscosity,
heat capacity, and conductivity give a wall-property correction of one;
temperature-dependent properties are needed to represent its variation.

This compares **coolant performance** with fixed gas-side conditions. Changing
the coolant name does not recalculate combustion chemistry. For another fuel
burned in the engine, also supply its gas properties in `main.py`.

To graph the IPA, ethanol, and saved RPA hot-wall temperature distributions,
select both coolants, run `python main.py`, then `python plot_temperatures.py`.
The graph is saved to `outputs/temperature_distributions.png`, with the full
chamber/nozzle and a throat-to-exit detail. RPA is labeled as the ethanol proxy.

## Model

The unchanged `contour.dxf` supplies the chamber/nozzle radius. Channel widths
vary with radius; velocity and hydraulic diameter are recalculated locally.
The gas side uses the textbook Bartz relation and constant-gamma area-Mach
flow, with recovery factor 0.90. Gas convection, wall conduction, and coolant
convection are solved together. Coolant flows from exit to injector.

The coolant side uses only **corrected Gnielinski**, with smooth-channel Darcy
friction and `(Pr_bulk / Pr_wall)^0.11`. Channel floors and sidewall fins are
included, normalized to hot-wall area. Bulk and coolant-wall temperatures are
iterated with the local properties. No other coolant correlation is included.

IPA density and viscosity use
[Kerscher et al., Eqs. 3 and 5](https://doi.org/10.1007/s10765-023-03294-z),
within their common 278-403 K range, without pressure corrections. Ethanol
density and viscosity use the installed RPA 2.3 liquid table interpolated to
3.447 MPa; its inlet heat capacity and conductivity are held constant. The
ethanol table is embedded in `Coolant.py`, so no separate property files are needed.

`wall_property_limit` caps the temperature used to evaluate wall properties;
it never clips the solved wall temperature. The retained IPA trial uses 403 K,
and the ethanol case uses 413 K. Set it to `None` to reject wall temperatures
outside the property range. Bulk and outlet properties always reject extrapolation.

Results are single-phase estimates. Near-transition flow, the wall-property
limit, constant heat capacity/conductivity, and unverified gas transport inputs
remain modeling limitations. Pressure loss, boiling, film cooling, radiation,
and axial conduction are omitted. Summary columns report limited wall-property
segments and segments below Re = 10,000.

Temperatures are K, coefficients W/(mÂ² K), heat flux W/mÂ², heat rate W, and
areas mÂ². CSV `x` and `radius_mm` are mm; coolant velocity is m/s.

## Files

| File | Purpose |
| --- | --- |
| `main.py` | Engine inputs, selected coolant runs, CSV results |
| `plot_temperatures.py` | Graph saved IPA, ethanol, and RPA wall temperatures |
| `Coolant.py` | Coolant inputs/properties, corrected Gnielinski, sidewall fins |
| `Contour.py` | DXF loading and contour geometry |
| `Geometry.py` | Cooling channels and heated segments |
| `Thermo.py` | Gas Mach number and temperatures |
| `Bartz.py` | Gas heat-transfer coefficient |
| `Wall.py` | Coupled wall balance and coolant heating |

`v1` preserves the original project and is not used by the model.
