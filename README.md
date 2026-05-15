# SAR Image Processing — Azimuth Compression

MATLAB implementation of azimuth compression for Synthetic Aperture Radar (SAR) data. The project takes raw, unfocused SAR signal data and produces a focused image through matched filtering in the azimuth direction, then refines it with multi-look processing.

## Background

A SAR sensor mounted on a moving platform (here, an aircraft) illuminates the ground with a wide beam. Each ground reflector is hit by many pulses as the aircraft passes by, producing a long, defocused response along the azimuth (flight) direction. Azimuth compression collapses that response back to a sharp point by convolving the raw data with a matched filter — the conjugate of the expected azimuth chirp. The result is a high-resolution image of the imaged scene.

## What the code does

- **Geometry & resolution analysis** — computes the synthetic aperture length, the required filter length, and theoretical azimuth resolution
- **Spectrum analysis** — estimates the azimuth bandwidth from the data and compares it to theory
- **Single-line focusing** — focuses one row of the raw signal to verify the filter design and tune the `fR/vx` parameter for best resolution
- **Full-image focusing** — applies a range-dependent matched filter bank to the full 350-row signal via frequency-domain convolution
- **Multi-look processing** — averages 9 looks to reduce speckle, then measures the resulting resolution

## File overview

| File | Purpose |
|------|---------|
| `Ans.m` | Main script — runs all tasks end to end |
| `sarfilter.m` | Applies the matched-filter bank to the signal (FFT-based, row by row) |
| `createHfilter.m` | Builds the range-dependent azimuth matched filter |
| `fftrows.m` / `ifftrows.m` | Row-wise FFT / IFFT helpers (zero-padded to 4096) |
| `spectrumplot.m` | Plots the averaged azimuth power spectrum |
| `filterplot.m` | Plots filter magnitude, phase, and group delay |
| `imageplot.m` | Displays SAR images with a fourth-root stretch for visibility |

## Data files

The raw signal (`sarlab.mat`) and reference image (`sarimage.mat`) are **not** included in this repository. They are course material and not mine to redistribute. To run the code, place those files in the project root.

## Running it

Open MATLAB, `cd` into the project folder, and run:

```matlab
Ans
```

The script will load the data, work through each task, and produce figures along the way.

## Key parameters

| Symbol | Value | Meaning |
|--------|-------|---------|
| ρ_min | 10446 m | Near-range distance |
| Δρ | 4 m | Range pixel size |
| λ_c | 0.0566 m | Carrier wavelength (C-band) |
| v_x | 131 m/s | Platform velocity |
| D | 1.3 m | Antenna length |
| Δx | 0.43 m | Azimuth sample spacing |
| f_R/v_x | ~2.33 m⁻¹ | PRF / velocity (tuned for best focus) |

## Notes

This was completed as a signal-processing course exercise. Posted here as a learning reference.
