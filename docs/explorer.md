---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Diurnal Cycle Explorer

Select an altitude, activate a chemical-family preset, then click individual
species in the legend to add or remove them. Grey bands mark nighttime (J = 0).
The horizontal axis follows DIURN's continuous noon-to-noon integration; the
clock labels therefore begin and end at 12:00.

```{code-cell} ipython3
:tags: [hide-input]
:execution_timeout: 600

from pratmo import Box, ChemistryOptions, IMPLICIT_SPECIES_NAMES, Model
import json
from IPython.display import HTML

# ── Species metadata ──────────────────────────────────────────────────────────

ALL_SPECIES = list(IMPLICIT_SPECIES_NAMES)

LABELS = {
    "no": "NO", "no2": "NO₂", "no3": "NO₃", "n2o5": "N₂O₅",
    "hno3": "HNO₃", "h": "H", "oh": "OH", "ho2": "HO₂",
    "h2o2": "H₂O₂", "o": "O(³P)", "o3": "O₃",
    "bro": "BrO", "br": "Br", "hbr": "HBr", "hno2": "HNO₂",
    "hcl": "HCl", "cl": "Cl", "cl2": "Cl₂", "clo": "ClO",
    "clono2": "ClONO₂", "hno4": "HNO₄", "hocl": "HOCl",
    "brono2": "BrONO₂", "hobr": "HOBr", "h2co": "H₂CO",
    "ch3o2": "CH₃O₂", "ch3o2h": "CH₃O₂H",
    "oclo": "OClO", "cl2o2": "Cl₂O₂", "brcl": "BrCl",
    # Iodine
    "i": "I", "io": "IO", "hoi": "HOI", "iono2": "IONO₂", "hi": "HI",
    "oio": "OIO", "i2": "I₂", "i2o2": "I₂O₂", "i2o3": "I₂O₃",
    "i2o4": "I₂O₄",
}

COLORS = {
    # NOx — reds/oranges
    "no":    "#e74c3c", "no2":   "#e67e22", "no3":  "#f39c12",
    "n2o5":  "#d35400", "hno3":  "#c0392b", "hno2": "#fa8072", "hno4": "#f0b27a",
    # HOx — blues
    "h":     "#85c1e9", "oh":    "#1a5276", "ho2":  "#2980b9", "h2o2": "#5dade2",
    # ClOx — greens
    "cl":    "#145a32", "clo":   "#1e8449", "cl2":  "#27ae60", "cl2o2": "#52be80",
    "hocl":  "#0e6655", "hcl":   "#0b5345", "clono2": "#82e0aa", "oclo": "#a9dfbf",
    # BrOx — purples
    "bro":   "#6c3483", "br":    "#7d3c98", "hbr":  "#8e44ad", "hobr": "#bb8fce",
    "brono2":"#5b2c6f", "brcl":  "#d2b4de",
    # Oxygen — browns
    "o":     "#784212", "o3":    "#9c640c",
    # Organics — pinks
    "h2co":  "#c71585", "ch3o2": "#e91e8c", "ch3o2h": "#ff69b4",
    # IOx — teals
    "i":     "#0077b6", "io":    "#00b4d8", "hoi":  "#48cae4",
    "iono2": "#023e8a", "hi":    "#90e0ef",
    "oio":   "#0096c7", "i2":    "#006d77", "i2o2": "#83c5be",
    "i2o3":  "#2a9d8f", "i2o4": "#264653",
}

PRESETS = {
    "NOx":      ["no", "no2", "no3", "n2o5", "hno3", "hno2", "hno4"],
    "HOx":      ["oh", "ho2", "h2o2", "h"],
    "ClOx":     ["cl", "clo", "cl2", "cl2o2", "hocl", "hcl", "clono2", "oclo"],
    "BrOx":     ["bro", "br", "hbr", "hobr", "brono2", "brcl"],
    "IOx":      ["i", "io", "hoi", "iono2", "hi", "oio", "i2", "i2o2", "i2o3", "i2o4"],
    "Ozone":    ["o", "o3"],
    "Organics": ["h2co", "ch3o2", "ch3o2h"],
}

DEFAULT_PRESET = "NOx"

# ── Run model for multiple altitudes ──────────────────────────────────────────

model = Model()
ALT_LEVELS = [8, 12, 18, 22]   # ~14, ~22, ~34, ~42 km
O3_SCALE_FACTOR = 1.0 # multiply initial box O3; lower this to test less ozone
NOY_SCALE_FACTOR = 1.0 # multiply initial total NOy; lower this to reduce NO2/IONO2
AEROSOL_AREA_SCALE_FACTOR = 1.0 # scales generic aerosol area for heterogeneous chemistry
AEROSOL_AREA_BY_LEVEL = {
    8: 1.0,
    12: 0.05,
    18: 0.0,
    22: 0.0,
}


payload = {}
for level in ALT_LEVELS:
    aerosol_area = AEROSOL_AREA_BY_LEVEL.get(level, 0.0) * AEROSOL_AREA_SCALE_FACTOR
    out = model.diurnal(
        latitude=0.0, day=120,
        boxes=[Box.at_level(
            level,
            aerosol_surface_area_um2_cm3=aerosol_area,
        )],
        chemistry=ChemistryOptions(iodine=True),
    )

    if O3_SCALE_FACTOR != 1.0 or NOY_SCALE_FACTOR != 1.0:
        init = out.boxes[0].long_lived
        init.o3 *= O3_SCALE_FACTOR
        init.noy *= NOY_SCALE_FACTOR
        out = model.diurnal(
            latitude=0.0, day=120,
            boxes=[Box.at_level(
                level,
                aerosol_surface_area_um2_cm3=aerosol_area,
            )],
            chemistry=ChemistryOptions(iodine=True),
            initial_mixing_ratios=[init],
        )

    ts  = out.time_series[0]
    snap = out.boxes[0]

    # DIURN is integrated continuously from noon to the following noon. Keep
    # that model order: sorting the cyclic HHMM labels splices the two distinct
    # noon boundary states together and creates an artificial discontinuity.
    steps_ordered = ts.steps
    hours = [step.elapsed_seconds / 3600.0 for step in steps_ordered]
    assert hours[0] == 0.0 and abs(hours[-1] - 24.0) < 1e-12
    assert all(a < b for a, b in zip(hours, hours[1:]))
    sp_data = {sp: [getattr(step.implicit, sp) for step in steps_ordered]
               for sp in ALL_SPECIES}

    # Fixed ticks every 2 h; labels show clock time for the noon-to-noon orbit.
    tick_vals = list(range(0, 25, 2))
    tick_text = [f"{(12 + h) % 24:02d}:00" for h in tick_vals]

    # Nighttime spans in decimal hours: contiguous runs where OH < 0.1 % of daily max
    oh = sp_data["oh"]
    threshold = max(oh) * 0.001
    night_spans, in_night, span_start_h = [], False, 0.0
    dt = [(hours[i] - hours[i-1]) / 2 for i in range(1, len(hours))]
    for i, v in enumerate(oh):
        half = dt[i-1] if i > 0 else 0.5
        if v < threshold and not in_night:
            span_start_h = hours[i] - half
            in_night = True
        elif v >= threshold and in_night:
            night_spans.append([span_start_h, hours[i-1] + half])
            in_night = False
    if in_night:
        # Extend the trailing night band to midnight so both nights are symmetric
        night_spans.append([span_start_h, 24.5])

    payload[str(level)] = {
        "alt_km":      round(snap.altitude_km, 0),
        "x":           hours,
        "tick_idx":    tick_vals,
        "tick_text":   tick_text,
        "species":     sp_data,
        "night_spans": night_spans,
        "o3_scale":    O3_SCALE_FACTOR,
    }

# ── Build widget HTML ─────────────────────────────────────────────────────────

alt_options = "\n".join(
    f'<option value="{lvl}">{int(payload[str(lvl)]["alt_km"])} km  (level {lvl})</option>'
    for lvl in ALT_LEVELS
)

preset_btns = "\n".join(
    f'<button type="button" class="px-btn" data-preset="{name}" aria-pressed="false" '
    f'onclick="dxSetPreset(\'{name}\')">{name}</button>'
    for name in PRESETS
)

widget = f"""
<style>
#dx-widget {{
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  margin: 8px 0 24px;
}}
#dx-widget .dx-bar {{
  display: flex;
  flex-wrap: wrap;
  gap: 12px;
  align-items: center;
  padding: 10px 14px;
  background: #f8f9fa;
  border: 1px solid #dee2e6;
  border-radius: 8px;
  margin-bottom: 10px;
  color: #212529;
}}
#dx-widget .dx-bar label {{
  font-weight: 600;
  font-size: 13px;
  color: #495057;
  margin-right: 4px;
  white-space: nowrap;
}}
#dx-widget select {{
  padding: 5px 8px;
  border: 1px solid #ced4da;
  border-radius: 5px;
  font-size: 13px;
  cursor: pointer;
  background: white;
  color: #212529;
}}
.dx-group {{
  display: flex;
  flex-wrap: wrap;
  gap: 5px;
  align-items: center;
}}
.px-btn {{
  padding: 4px 11px;
  border: 1px solid #adb5bd;
  border-radius: 16px;
  background: white;
  font-size: 12px;
  cursor: pointer;
  color: #495057;
  transition: background 0.12s, color 0.12s, border-color 0.12s;
  white-space: nowrap;
}}
.px-btn:hover {{ background: #e9ecef; border-color: #6c757d; }}
.px-btn.active {{
  background: #0d6efd;
  border-color: #0d6efd;
  color: #fff;
  font-weight: 600;
}}
.px-btn.extra {{ border-style: dashed; }}
</style>

<div id="dx-widget">
  <div class="dx-bar">
    <div>
      <label for="dx-alt">Altitude:</label>
      <select id="dx-alt" onchange="dxUpdateAlt(this.value)">{alt_options}</select>
    </div>
    <div>
      <label for="dx-scale">Scale:</label>
      <select id="dx-scale" onchange="dxUpdateScale(this.value)">
        <option value="log" selected>Log</option>
        <option value="linear">Linear</option>
      </select>
    </div>
    <div class="dx-group">
      <label>Preset:</label>
      {preset_btns}
      <button type="button" class="px-btn extra" data-preset="_all" aria-pressed="false" onclick="dxSetPreset('_all')">All</button>
      <button type="button" class="px-btn extra" data-preset="_none" aria-pressed="false" onclick="dxSetPreset('_none')">None</button>
    </div>
  </div>
  <div id="dx-chart"></div>
</div>

<script src="https://cdn.plot.ly/plotly-2.26.0.min.js" charset="utf-8"></script>
<script>
(function () {{
  var DATA     = {json.dumps(payload)};
  var PRESETS  = {json.dumps(PRESETS)};
  var COLORS   = {json.dumps(COLORS)};
  var LABELS   = {json.dumps(LABELS)};
  var SPECIES  = {json.dumps(ALL_SPECIES)};
  var DEFAULT  = "{DEFAULT_PRESET}";

  var curAlt    = "{ALT_LEVELS[0]}";
  var curPreset = DEFAULT;
  var curScale  = "log";

  // ── helpers ──────────────────────────────────────────────────────────────

  function activeSet(preset) {{
    if (preset === "_all")  return new Set(SPECIES);
    if (preset === "_none") return new Set();
    return new Set(PRESETS[preset] || []);
  }}

  function nightShapes(altKey) {{
    return DATA[altKey].night_spans.map(function (s) {{
      return {{
        type: "rect", xref: "x", yref: "paper",
        x0: s[0], x1: s[1], y0: 0, y1: 1,
        fillcolor: "#c8c8c8", opacity: 0.30,
        line: {{ width: 0 }}, layer: "below",
      }};
    }});
  }}

  function makeTraces(altKey, preset) {{
    var d   = DATA[altKey];
    var act = activeSet(preset);
    return SPECIES.map(function (sp) {{
      return {{
        x: d.x,
        y: d.species[sp],
        name: LABELS[sp],
        type: "scatter",
        mode: "lines",
        line: {{ color: COLORS[sp], width: 2.5 }},
        visible: act.has(sp) ? true : "legendonly",
        hovertemplate: "<b>" + LABELS[sp] + "</b>: %{{y:.2e}} cm⁻³<extra></extra>",
      }};
    }});
  }}

  function makeLayout(altKey) {{
    var d = DATA[altKey];
    return {{
      margin: {{ t: 10, r: 200, b: 75, l: 90 }},
      xaxis: {{
        title: {{ text: "Local solar time (noon-to-noon)", standoff: 10 }},
        tickmode: "array",
        tickvals: d.tick_idx,
        ticktext: d.tick_text,
        tickangle: -40,
        tickfont: {{ size: 11 }},
        showgrid: true, gridcolor: "#f0f0f0",
        range: [-0.5, 24.5],
      }},
      yaxis: {{
        title: {{ text: "Number density (cm⁻³)", standoff: 10 }},
        type: curScale,
        exponentformat: "power",
        showgrid: true, gridcolor: "#f0f0f0",
      }},
      legend: {{
        x: 1.02, y: 1, xanchor: "left", yanchor: "top",
        font: {{ size: 11 }},
        bgcolor: "rgba(255,255,255,0.9)",
        bordercolor: "#dee2e6", borderwidth: 1,
      }},
      hovermode: "x",
      height: 520,
      plot_bgcolor: "#ffffff",
      paper_bgcolor: "#ffffff",
      shapes: nightShapes(altKey),
    }};
  }}

  // ── public API ────────────────────────────────────────────────────────────

  window.dxUpdateAlt = function (altKey) {{
    curAlt = altKey;
    var d   = DATA[altKey];
    var upd = {{
      x: SPECIES.map(function () {{ return d.x; }}),
      y: SPECIES.map(function (sp) {{ return d.species[sp]; }}),
    }};
    Plotly.restyle("dx-chart", upd);
    Plotly.relayout("dx-chart", {{
      "xaxis.tickvals": d.tick_idx,
      "xaxis.ticktext": d.tick_text,
      "xaxis.range": [-0.5, 24.5],
      shapes: nightShapes(altKey),
    }});
  }};

  window.dxSetPreset = function (preset) {{
    curPreset = preset;
    var act     = activeSet(preset);
    var visible = SPECIES.map(function (sp) {{ return act.has(sp) ? true : "legendonly"; }});
    Plotly.restyle("dx-chart", {{ visible: visible }});
    document.querySelectorAll(".px-btn").forEach(function (b) {{
      var active = b.dataset.preset === preset;
      b.classList.toggle("active", active);
      b.setAttribute("aria-pressed", active ? "true" : "false");
    }});
  }};

  window.dxUpdateScale = function (scale) {{
    curScale = scale;
    Plotly.relayout("dx-chart", {{
      "yaxis.type": curScale,
    }});
  }};

  // ── init ──────────────────────────────────────────────────────────────────

  function init() {{
    Plotly.newPlot(
      "dx-chart",
      makeTraces(curAlt, curPreset),
      makeLayout(curAlt),
      {{ responsive: true, displaylogo: false,
         modeBarButtonsToRemove: ["lasso2d", "select2d"] }}
    );
    dxSetPreset(DEFAULT);
    dxUpdateAlt(curAlt);
  }}

  if (typeof Plotly !== "undefined") {{
    init();
  }} else {{
    var s = document.querySelector('script[src*="plotly"]');
    if (s) s.addEventListener("load", init);
  }}
}})();
</script>
"""

HTML(widget)
```
