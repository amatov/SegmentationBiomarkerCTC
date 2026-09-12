## Quick start

This repository covers segmentation of PSMA-labeled circulating tumor
cells (CTCs) at two resolutions, plus a set of general-purpose supporting
tools. See [DEPENDENCIES.md](DEPENDENCIES.md) for the MATLAB Compiler
Runtime requirement of the compiled executables.

## Repository contents

- **`Coverslip/`** -- the main analysis code, covering both coverslip and
  Gedi-chip samples:
  - `high_resolution/` -- `WaveletSeedingWatershedActiveContour/` (the
    wavelet-seeding + active-contour + watershed segmentation algorithm).
  - `low_resolution/` -- `find_ctc.exe` / `optimal_ctc.exe` (compiled
    low-resolution whole-slide analysis), `codeCTC.zip`, `genOptimalCTC2.m`,
    the `.ini` configuration files (covering both coverslip and Gedi-chip
    datasets), and `run_optimal_detection_9283.bat`.
  - `Segmentation3D/` -- 3D segmentation code.
  - `spotDetection/` -- bright-spot detection via stationary wavelet
    transform.
  - The Matlab scripts at this folder's top level (`CTCQ_*.m`,
    `Coverslips_new_code.m`) are the high-resolution cropped-image
    analysis.
- **`additional_tools/`** -- general-purpose tools used across the
  project, not specific to coverslip/chip samples:
  - `SVM/` -- Support Vector Machine classification code.
  - `HoughTransform/` -- Hough Transform line/circle detection code (see
    that folder's own README for its specific license terms).
- **`media/`** -- supplementary video, image, and PDF files.
- **License:** see [LICENSE](LICENSE) -- research/educational use, with
  separate terms noted for bundled third-party components.

## About

Algorithm for the segmentation of prostate-specific membrane antigen labeling in ciruculating tumor cells of metastatic patients I designed and partially implemented

The patient blood samples analyzed were from clinical studies with IRB protocols 0804009740 and 0707009283 at Cornell Medicine and NCT01718353 phase II clinical trial sponsored by Sanofi for early switch from Docetaxel to Cabazitaxel during the treatment of metastatic castrate-resistant prostate cancer; see the examples of image processing of whole slides in folder Coverslip

Watch my CTC presentation, the second part of this seminar, at the University of Central Florida in 2013: https://youtu.be/kTYyltX9RFg?t=1564

For detailed information, see: https://www.researchgate.net/publication/387438061_Microtubule_Regulation_in_Cancer_Cells
