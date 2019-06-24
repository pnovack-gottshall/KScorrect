# KScorrect 1.4.0 (2019-06-30)
* Added new function 'ks_test_stat' to more quickly calculate test statistic D.
* Added new argument 'varModel' to 'LcKS' to allow users to specify 'equal' or 'variable' variance models for mixture distributions.
* Added option to run 'LcKS' in parallel, using 'doParallel' and 'foreach' infrastructure. Added new parameters 'parallel' and 'cores' to control parallel behavior.

# KScorrect 1.3.0 (2019-01-05)
* Changed spline method from default to 'hyman' in qmixnorm that is more appropriate for handling monotonic sequences, and improves performance under certain mixture models. (Thanks to Qiong Zhang for catching the error.)

# KScorrect 1.2.5 (2018-09-17)
* Changed number of randomly generated samples per mixture components ('nr') to 1000, to match help file text. nr = 1000 is much faster and nearly equally precise.

# KScorrect 1.2.4 (2018-08-14)
* Silenced 'verbose' arg when calling 'mclust' functions.

# KScorrect 1.2.3 (2018-02-23)
* Reformatted functions to Google / Hadley style.

# KScorrect 1.2.2 (2017-04-23)
* Clarify text so better explains behavior based on 'mclust' package functions.

# KScorrect 1.2.1 (2017-01-10)
* Improved text, with corrected textual errors.

# KScorrect 1.2.0 (2016-03-19)
* Text clarifications to Help files prior to posting to CRAN.
* Removed calculation of standard error in LcKS.

# KScorrect 1.1.0
* Steve Wang re-wrote qmixnorm to improve performance and make the function more compatible with behavior in other quantile functions, such as qnorm.
* Improved handling of extreme probabilities in qmixnorm, based on algorithm by Luca Scrucca.
* Corrected behavior when single 'sd' specified in mixture model functions so that results in equal-variance mixture model, and calls warning to alert to behavior

# KScorrect 1.0.0
* First release (to GitHub).
