## Test Environments

* Local macOS 27.2 install, R 4.6.2
* win-builder, R-devel (2026-11-05 r90840 ucrt)

## R CMD check results

0 errors | 0 warnings | 0 notes

This holds on both test environments above.

## Update of existing CRAN package

* This is a patch on an existing CRAN package to address the "Lost braces" NOTE in the Rd files currently appearing in the CRAN checks.
* No changes to the package code.

## Reverse dependencies

* We checked the 2 reverse dependencies (KrigInv, pGPx) against this version, with `R CMD check --run-donttest`: no new problems. Both currently show a "Lost braces" Rd NOTE, which is unrelated to anMC and will be fixed in their upcoming patch releases.
