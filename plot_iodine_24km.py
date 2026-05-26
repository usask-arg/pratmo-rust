import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pratmo import PratmoModel, DiurnConfig, DiurnBoxSpec

m = PratmoModel.with_defaults()
out = m.run_diurn(DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,
    integration_days=20,
    boxes=[DiurnBoxSpec(altitude_level=13)],  # 24 km
))

dm = out.boxes[0].air_density_cm3  # air number density cm⁻³
alt_km = out.boxes[0].altitude_km
print(f"altitude = {alt_km:.1f} km,  dm = {dm:.3e} cm⁻³")

steps = out.time_series[0].steps
time_hr = []
io_ppt = []
hoi_ppt = []
iono2_ppt = []
i_ppt = []
hi_ppt = []

for s in steps:
    hhmm = s.time_hhmm
    hr = (hhmm // 100) + (hhmm % 100) / 60.0
    # wrap midnight so timeline is continuous (6 am to 6 am next day)
    if hr < 6:
        hr += 24
    time_hr.append(hr)
    imp = s.implicit
    ppt = 1e12 / dm
    io_ppt.append(imp.io * ppt)
    hoi_ppt.append(imp.hoi * ppt)
    iono2_ppt.append(imp.iono2 * ppt)
    i_ppt.append(imp.i * ppt)
    hi_ppt.append(imp.hi * ppt)

time_hr = np.array(time_hr)
io_ppt = np.array(io_ppt)
hoi_ppt = np.array(hoi_ppt)
iono2_ppt = np.array(iono2_ppt)
i_ppt = np.array(i_ppt)
hi_ppt = np.array(hi_ppt)

fig, ax = plt.subplots(figsize=(10, 5))

# night shading: before 6am and after 18hr (6pm)
ax.axvspan(6, 18, alpha=0.12, color="gold", label="Daytime")
ax.axvspan(6, 6, alpha=0)   # force x limits inclusion

ax.plot(time_hr, io_ppt,    color="royalblue",   lw=2,   label="IO")
ax.plot(time_hr, hoi_ppt,   color="darkorange",  lw=1.5, label="HOI")
ax.plot(time_hr, iono2_ppt, color="green",       lw=1.5, label="IONO₂")
ax.plot(time_hr, i_ppt,     color="purple",      lw=1.5, label="I")
ax.plot(time_hr, hi_ppt,    color="red",         lw=1.5, label="HI")

xticks = [6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
ax.set_xticks(xticks)
ax.set_xticklabels([f"{t%24:02d}:00" for t in xticks])
ax.set_xlim(6, 30)
ax.set_xlabel("Local time (UTC)", fontsize=12)
ax.set_ylabel("Mixing ratio (ppt)", fontsize=12)
ax.set_title(f"Iodine diurnal cycle at {alt_km:.0f} km  (0°N, Julian day 120)", fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig("/tmp/iodine_24km.png", dpi=150)
print("Saved /tmp/iodine_24km.png")
print(f"\nPeak IO  = {io_ppt.max():.4f} ppt")
print(f"Night IO = {io_ppt.min():.2e} ppt")
print(f"HOI range: {hoi_ppt.min():.4f} – {hoi_ppt.max():.4f} ppt")
print(f"IONO2 range: {iono2_ppt.min():.4f} – {iono2_ppt.max():.4f} ppt")
