## Resubmission

This is a resubmission of autoRasch, which was archived on 2025-11-27 because
issues were not corrected in time. The CRAN reminders were sent to the
previous maintainer address (f.wijayanto@cs.ru.nl), which is no longer active.
The maintainer (Feri Wijayanto) is unchanged; only the contact address has been
updated to feri.wijayanto@uii.ac.id.

The following issues have been fixed:

* Removed `CXX_STD = CXX11` (RcppArmadillo fallback compilation warning and
  C++ specification note).
* Added the `autoRasch-package` alias via the `"_PACKAGE"` roxygen2 sentinel.
* Fixed the Rd "Lost braces in \itemize" notes.
* Renamed the data files with the `.Rdata` extension to `.RData`.

## Test environments

* local macOS (aarch64), R 4.6.1
* TODO before submission: win-builder (R-devel), rhub (Linux R-devel, Windows)

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: New submission / package was archived on CRAN (see above).
