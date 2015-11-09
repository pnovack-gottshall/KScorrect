## KScorrect
R package for Lilliefors-corrected Kolmogorov-Smirnoff goodness-of-fit tests

KScorrect implements the Lilliefors-corrected Kolmogorov-Smirnoff test for use in goodness-of-fit tests, suitable when using sample statistics as estimates of population parameters. It uses a resampling algorithm to estimate p-values from the resampling distribution. Coded to complement stats::ks.test, it is able to be used with a variety of continuous distributions, including normal, lognormal, univariate mixtures of normals, uniform, loguniform, exponential, gamma, and Weibull distributions.

Functions are also provided to generate random numbers and calculate density, distribution, and quantile functions are provided for use with the loguniform and mixture distributions.

## References

The most recent commit is currently: [![Travis-CI Build Status](https://travis-ci.org/pnovack-gottshall/KScorrect.svg?branch=master)](https://travis-ci.org/pnovack-gottshall/KScorrect) (Travis CI)

This package is authored by Phil Novack-Gottshall (<mailto:pnovack-gottshall@ben.edu>) and Steve C. Wang (<mailto:stevecwang@swarthmore.edu>) and offered under CC0.
