# A549 hpDPC worked example

This directory contains a compact example for reproducing the hpDPC spectral-fusion step from precomputed pDPC and TIE phase maps.

The example is intentionally limited to one representative sample. Additional experimental examples and the classification dataset will be archived separately in the associated Zenodo dataset.

## Files

```text
example/
├── README.md
├── run_example.m
└── A549/
    ├── A549_pDPC.tif
    ├── A549_TIE.tif
    └── A549_hpDPC_reference.tif
```

### `A549_pDPC.tif`

Precomputed pDPC phase reconstruction used as the higher-spatial-frequency fusion input.

### `A549_TIE.tif`

Precomputed TIE phase reconstruction used as the low-spatial-frequency fusion input.

### `A549_hpDPC_reference.tif`

Reference hpDPC reconstruction generated using the manuscript fusion parameters and stored as a 32-bit floating-point TIFF.

### `run_example.m`

Loads the supplied phase maps, calls `fusion/hpDPC_fusion.m`, and compares the recreated hpDPC phase map with the stored reference.

## Input data

The three supplied phase maps are:

- single-channel;
- 32-bit IEEE floating point;
- 2048 × 2448 pixels;
- expressed in radians.

The pDPC and TIE inputs are already reconstructed and co-registered on the same computational grid. The example does not perform pDPC reconstruction, TIE reconstruction, registration, normalization, or resizing.

## Fusion parameters

The example uses:

```text
dx      = 3.45e-6 / 20 m/pixel
k0      = 0.05
alpha   = 2
reg_eps = 1e-3
```

These values are defined near the top of `run_example.m` and can be modified for exploratory use.

## Running the example

From MATLAB:

```matlab
cd example
run_example
```

The script automatically adds the repository `fusion/` directory to the MATLAB path.

It reports:

- exact equality after conversion to the stored single-precision representation;
- root-mean-square error (RMSE);
- mean absolute error (MAE);
- maximum absolute error;
- RMSE normalized by the phase range of the reference image.

Small nonzero differences can occur across MATLAB versions or platforms because the fusion uses FFT-based floating-point calculations. The numerical error metrics should therefore be used in addition to the exact-match flag.

## Data license

The TIFF files in `A549/` are released under **CC BY 4.0**. See `../DATA_LICENSE.md`.

## Full supporting dataset

The larger Zenodo dataset will contain additional experimental examples, classification source data, predictions, summary metrics, and metadata.

**Dataset DOI:** to be added when the dataset is released.
