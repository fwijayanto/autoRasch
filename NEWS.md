# autoRasch 0.2.4

* Includes the changes of the unreleased development version 0.2.3
  (new datasets `DIFdataset` and `shortDSF`).
* Updated the maintainer's contact email address.
* Removed the `CXX_STD = CXX11` specification from `src/Makevars` and
  `src/Makevars.win`, which caused the RcppArmadillo "fallback compilation"
  warning and the "C++ specification" note.
* Fixed the Rd `\value` sections of `fitStats()`, `checkRel()` and the
  `pcm`/`pcm_dif` documentation ("Lost braces" notes) and the `\describe`
  list in `autoRaschOptions()`.
* Documented the package with the `"_PACKAGE"` sentinel so that the
  `autoRasch-package` alias is available (`?"autoRasch-package"`).
* Renamed the data files `DIFdataset.Rdata` and `shortDSF.Rdata` to use the
  `.RData` extension so that they are recognised by R.
* Added `URL` and `BugReports` fields; updated the vignette author.

# autoRasch 0.2.2

* See the previous release on CRAN.
