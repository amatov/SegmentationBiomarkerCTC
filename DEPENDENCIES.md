# Dependencies

## `coverslip/low_resolution/find_ctc.exe`

Requires the **Matlab Compiler Runtime (MCR) version 7.16**. Neither the
MCR nor `MCRInstaller.exe` is included in this repository -- see
`coverslip/low_resolution/README.md` for the original deployment notes.
Download the matching MCR installer from MathWorks before running
`find_ctc.exe`.

## `coverslip/low_resolution/optimal_ctc.exe` / `awt.mexw64`

Also Matlab-Compiler-based binaries; same MCR requirement applies.
`awt.mexw64` is a compiled MEX file built from `awt.cpp` -- to rebuild it
from source on a different platform, use Matlab's `mex` compiler.

## `coverslip/low_resolution/codeCTC.zip`

Packaged separately; see its own contents for build/runtime requirements
once unzipped.

## Matlab scripts (`.m` files throughout the repository)

No specific toolbox versions are recorded. Scripts using `activecontour`,
`watershed`, `regionprops`, and similar functions require Matlab's Image
Processing Toolbox.
