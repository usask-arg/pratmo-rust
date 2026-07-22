"""Interactive Plotly helpers for PRATMO outputs and atmospheric inputs."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from pratmo.quantities import mixing_ratio_as, number_density_as


def _plotly():
    try:
        import plotly.graph_objects as go
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "Interactive plotting requires Plotly. Install it with `pip install 'pratmo[plot]'`."
        ) from exc
    return go


def plot_profile(
    output,
    species: str | Sequence[str] = ("o3", "no2", "oh"),
    *,
    kind: str = "number_density",
    unit: str | None = None,
    log_x: bool = True,
):
    """Plot one or more vertical chemistry, mixing-ratio, or J-value profiles."""

    go = _plotly()
    names = [species] if isinstance(species, str) else list(species)
    key = kind.strip().lower().replace("-", "_")
    if key in {"number_density", "implicit", "species"}:
        requested_unit = unit or "cm-3"
        getter = output.species_profile
        convert = lambda value: number_density_as(value, requested_unit)
        label = f"Number density ({requested_unit})"
    elif key in {"mixing_ratio", "long_lived", "vmr"}:
        requested_unit = unit or "ppbv"
        getter = output.long_lived_profile
        convert = lambda value: mixing_ratio_as(value, requested_unit)
        label = f"Mixing ratio ({requested_unit})"
    elif key in {"jvalue", "j_value", "photolysis"}:
        if unit not in {None, "s-1", "s^-1"}:
            raise ValueError("J-values are available in s-1")
        getter = output.jvalue_profile
        convert = np.asarray
        label = "Photolysis rate (s⁻¹)"
    else:
        raise ValueError("kind must be number_density, mixing_ratio, or jvalue")

    figure = go.Figure()
    for name in names:
        figure.add_scatter(
            x=convert(getter(name)),
            y=output.altitude_km,
            mode="lines+markers",
            name=name.upper(),
            hovertemplate=f"{name.upper()}: %{{x:.3g}}<br>Altitude: %{{y:.2f}} km<extra></extra>",
        )
    figure.update_layout(
        xaxis_title=label,
        yaxis_title="Altitude (km)",
        template="plotly_white",
        hovermode="closest",
        legend_title_text="Field",
    )
    figure.update_xaxes(type="log" if log_x else "linear")
    return figure


def plot_diurnal(
    output,
    species: str | Sequence[str] = ("oh", "ho2", "no2"),
    *,
    box: int = 0,
    unit: str = "cm-3",
    log_y: bool = True,
):
    """Plot resolved species through the noon-to-noon DIURN trajectory."""

    go = _plotly()
    names = [species] if isinstance(species, str) else list(species)
    if not 0 <= box < len(output):
        raise IndexError(f"box must be between 0 and {len(output) - 1}")
    hours = output.elapsed_seconds / 3600.0
    figure = go.Figure()
    for name in names:
        values = number_density_as(output.species_grid(name)[box], unit)
        figure.add_scatter(
            x=hours,
            y=values,
            mode="lines+markers",
            name=name.upper(),
            hovertemplate=f"{name.upper()}: %{{y:.3g}} {unit}<br>Elapsed: %{{x:.2f}} h<extra></extra>",
        )
    tick_values = np.arange(0.0, 25.0, 3.0)
    figure.update_layout(
        xaxis_title="Local solar time (noon to noon)",
        yaxis_title=f"Number density ({unit})",
        template="plotly_white",
        hovermode="x unified",
        legend_title_text=f"{output.altitude_km[box]:.1f} km",
    )
    figure.update_xaxes(
        tickmode="array",
        tickvals=tick_values,
        ticktext=[f"{int((12 + hour) % 24):02d}:00" for hour in tick_values],
    )
    figure.update_yaxes(type="log" if log_y else "linear")
    return figure


def plot_atmosphere(atmosphere):
    """Plot the pressure, temperature, ozone, and optional aerosol profile."""

    go = _plotly()
    from plotly.subplots import make_subplots

    if atmosphere.altitude_km is None:
        vertical = np.arange(len(atmosphere))
        vertical_label = "Profile level"
    else:
        vertical = atmosphere.altitude_km
        vertical_label = "Altitude (km)"
    columns = 4 if atmosphere.aerosol_surface_area_um2_cm3 is not None else 3
    figure = make_subplots(rows=1, cols=columns, shared_yaxes=True)
    fields = [
        (atmosphere.temperature_k, "Temperature (K)", False),
        (atmosphere.pressure_mb, "Pressure (hPa)", True),
        (
            mixing_ratio_as(atmosphere.ozone, "ppmv")
            if atmosphere.ozone_kind == "mixing_ratio"
            else atmosphere.ozone,
            "Ozone (ppmv)" if atmosphere.ozone_kind == "mixing_ratio" else "Ozone (cm⁻³)",
            True,
        ),
    ]
    if atmosphere.aerosol_surface_area_um2_cm3 is not None:
        fields.append(
            (
                atmosphere.aerosol_surface_area_um2_cm3,
                "Aerosol (µm² cm⁻³)",
                True,
            )
        )
    for column, (values, label, use_log) in enumerate(fields, start=1):
        figure.add_scatter(
            x=values,
            y=vertical,
            mode="lines+markers",
            showlegend=False,
            row=1,
            col=column,
        )
        figure.update_xaxes(title_text=label, type="log" if use_log else "linear", row=1, col=column)
    figure.update_yaxes(title_text=vertical_label, row=1, col=1)
    figure.update_layout(template="plotly_white", hovermode="closest")
    return figure


__all__ = ["plot_atmosphere", "plot_diurnal", "plot_profile"]
