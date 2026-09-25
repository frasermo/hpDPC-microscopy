# hpDPC microscopy

MATLAB code and a representative worked example accompanying the manuscript:

**Hybrid Transport-of-Intensity and Polarization Differential Phase Contrast for Extended Spatial-Frequency Phase Imaging**

This repository contains the hpDPC spectral-fusion implementation, the simulation used to illustrate the method, and one compact A549 example for reproducing the fusion step from precomputed phase maps.

## Scope

The public fusion code operates on **precomputed pDPC and TIE phase reconstructions**. 

The repository provides:

- the hpDPC Fourier-domain spectral-fusion function;
- a self-contained MATLAB simulation of the hpDPC concept;
- one representative A549 example with pDPC and TIE phase inputs;
- a stored hpDPC reference reconstruction for numerical comparison.

A separate Zenodo dataset will contain the larger supporting data package, including additional experimental examples, classification source data, per-patch predictions, summary metrics, and associated metadata. Those larger data are intentionally not duplicated in this GitHub repository.

## Repository structure

```text
hpDPC-microscopy/
├── README.md
├── LICENSE
├── DATA_LICENSE.md
├── CITATION.cff
├── .gitignore
│
├── fusion/
│   └── hpDPC_fusion.m
│
├── simulations/
│   └── hpdpc_simulation.m
│
└── example/
    ├── README.md
    ├── run_example.m
    └── A549/
        ├── A549_pDPC.tif
        ├── A549_TIE.tif
        └── A549_hpDPC_reference.tif
```

## Software requirements

The code was tested with **MATLAB R2024b**.

The fusion function, worked example, and simulation use MATLAB functionality only.

## hpDPC spectral fusion

The fusion implementation is provided in:

```text
fusion/hpDPC_fusion.m
```

The function combines precomputed TIE and pDPC phase maps in the Fourier domain using complementary radial weights. In the manuscript implementation, TIE supplies the low-spatial-frequency contribution and pDPC supplies the complementary higher-spatial-frequency contribution.

The two input phase maps must have identical dimensions and should be co-registered, expressed with a consistent phase sign convention, and sampled on the same computational grid before fusion.

A minimal call is:

```matlab
phiHpDPC = hpDPC_fusion( ...
    'phiTIE', phiTIE, ...
    'phiDPC', phiDPC, ...
    'dx', dx, ...
    'k0', 0.05, ...
    'alpha', 2, ...
    'reg_eps', 1e-3);
```

See the function header in `fusion/hpDPC_fusion.m` for the complete parameter description.

## Worked A549 example

The example is located in:

```text
example/
```

Run from MATLAB with:

```matlab
cd example
run_example
```

The script loads the supplied pDPC and TIE phase maps, reconstructs hpDPC using `fusion/hpDPC_fusion.m`, and compares the result with the stored reference reconstruction. It reports exact equality together with RMSE, mean absolute error, and maximum absolute error so that small floating-point differences can be identified explicitly.

The supplied TIFF files are:

- single-channel;
- 32-bit IEEE floating point;
- 2048 × 2448 pixels;
- expressed in radians.

See `example/README.md` for details.

## Simulation

The manuscript simulation is provided in:

```text
simulations/hpdpc_simulation.m
```

Run with:

```matlab
cd simulations
hpdpc_simulation
```

The script simulates and compares ground-truth phase, TIE, matched pDPC, unmatched pDPC, a hybrid TIE-pDPC reconstruction, and hpDPC. The partially coherent forward model uses discrete source-point sampling, with `targetSourceSamples = 80` for the manuscript configuration.

An optional source-sampling convergence test can be enabled inside the script.

## Supporting dataset

A separate Zenodo dataset will archive the larger supporting data package. It is expected to include:

- additional representative experimental phase data;
- pDPC, TIE, and hpDPC classification source data;
- per-patch labels and prediction probabilities;
- classification summary metrics;
- metadata and data dictionaries;
- additional publication examples.

**Dataset DOI:** to be added when the dataset is released.

## Citation

This repository accompanies the manuscript:

> F. Montandon, J. Knopp, C. A. Jacobs, and F. Nicolls, “Hybrid Transport-of-Intensity and Polarization Differential Phase Contrast for Extended Spatial-Frequency Phase Imaging,” submitted for publication, 2026.

A version-specific Zenodo DOI for the software release will be added when the release is archived.

The final journal citation and article DOI will be added after publication. Machine-readable citation metadata are provided in `CITATION.cff`.

## Versioning and archival release

The manuscript-associated software release should be tagged as:

```text
v1.0.0
```

The tagged release can then be archived on Zenodo. The version-specific Zenodo DOI should be used when citing the software version associated with the manuscript.

## Licensing

Original source code and software documentation in this repository are released under the **MIT License**. See `LICENSE`.

The example TIFF data are released under **Creative Commons Attribution 4.0 International (CC BY 4.0)**. See `DATA_LICENSE.md`.

## Authors

- Fraser Montandon †
- Jasmin Knopp
- Caron A. Jacobs
- Fred Nicolls

University of Cape Town, Cape Town, South Africa.

Additional affiliations are given in the associated manuscript.
† indicates the corresponding author.

## Contact

For questions about the software or supporting data, please contact the corresponding author or open an issue in this repository.
