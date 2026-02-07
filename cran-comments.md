## Resubmission

This is a resubmission. In this version I have:

* Saved and restored `par()` in all vignette plotting chunks via
  `oldpar <- par(...)` / `par(oldpar)`.
* Removed `<<-` usage from vignettes and tests; replaced with
  explicit environment-based trace collection.
* Removed unnecessary single quotes around technical terms in DESCRIPTION;
  added `inst/WORDLIST` for the spell checker instead.
* Spelled out "MLEs" as "maximum likelihood estimators" in DESCRIPTION.
* Removed all commented-out code from `@examples` sections (inline expected-value
  comments and standalone comment lines) across all R source files and regenerated
  man pages.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* Local: Ubuntu 24.04.3 LTS, R 4.3.3
* GitHub Actions: Ubuntu (release, devel), macOS (release), Windows (release)

## Downstream dependencies

This is a new package with no reverse dependencies.
