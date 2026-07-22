# Standard atmospheric levels

Standard levels are 1-based. Their pressure coordinate is fixed, while
temperature and hydrostatically derived altitude vary with the selected
latitude and season. The table below is a practical orientation for the
default case, 45°N on day 80; it is not a universal height lookup.

| Level | Approx. altitude (km) | Pressure (hPa) | Temperature (K) | Typical role |
|---:|---:|---:|---:|---|
| 1 | 0.0 | 1000 | 275.4 | Lower boundary |
| 5 | 8.6 | 316 | 228.9 | Upper troposphere |
| 8 | 14.1 | 133 | 217.3 | Lower stratosphere |
| 10 | 17.8 | 75.0 | 215.8 | Default CTM profile |
| 12 | 21.4 | 42.2 | 216.4 | Lower stratosphere |
| 15 | 27.0 | 17.8 | 221.8 | Default CTM and DIURN |
| 16 | 28.8 | 13.3 | 224.4 | Middle stratosphere |
| 18 | 32.7 | 7.50 | 232.5 | Middle stratosphere |
| 20 | 36.7 | 4.22 | 243.7 | Default CTM profile |
| 22 | 40.9 | 2.37 | 254.3 | Upper stratosphere |
| 24 | 45.2 | 1.33 | 261.8 | Upper stratosphere |
| 25 | 47.5 | 1.00 | 263.2 | Default CTM profile |
| 30 | 58.3 | 0.237 | 246.6 | Lower mesosphere |

```python
from pratmo import Box, Model

output = Model().ctm(
    boxes=[Box.at_level(level) for level in (8, 12, 16, 20, 24)]
)
print(output.altitude_km)  # use the run's coordinates, not the table above
print(output.pressure_mb)
```

The embedded grid extends to level 41. Most documented examples stay in the
stratospheric part of the grid where this model is intended to be used. A
custom atmosphere should be used when exact instrument altitudes or a specific
retrieval grid matter.
