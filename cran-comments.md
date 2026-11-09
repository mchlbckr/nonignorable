## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

  The note also lists "Debiased", "Nonignorable" and "cker" as possibly
  misspelled words in DESCRIPTION. The first two are technical terms;
  "cker" is part of the author's surname, Bücker.

## Test environments

* local: macOS, R release
* win-builder: R-devel, R release
* GitHub Actions: macOS and Windows (R release), Ubuntu (R oldrel-1,
  release and devel)

## Notes for the reviewer

* The vignette "applications" computes one example live from the suggested
  packages 'speff2trial' and 'mgcv'. It checks for both with
  requireNamespace() and skips these chunks when either is missing.
* Examples, tests and vignettes run in about 5, 45 and 16 seconds on the
  local machine.
