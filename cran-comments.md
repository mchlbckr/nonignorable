## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

## Test environments

* local: macOS, R release
* win-builder: R-devel

## Notes for the reviewer

* The vignette "applications" computes one example live from the suggested
  packages 'speff2trial' and 'mgcv'. It checks for both with
  requireNamespace() and skips these chunks when either is missing.
* Examples, tests and vignettes run in about 5, 45 and 16 seconds on the
  local machine.
