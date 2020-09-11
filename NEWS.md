# adopter 0.4.0

* included JSS article as vignette

# adoptr 0.3.2

* bugfix for binomial pdf/cdf

# adoptr 0.3.1

* added support for binomial endpoint
* initial design creation more convenient
* extended summary() function for two stage designs

# adoptr 0.3.0

* bugfix for pretty printing of optimization results
* streamlined code for internal evaluation of integral scores during optimization; 
substantial speedup!

# adoptr 0.2.3

* reworked the class printing system to be more informative
* minimize() now also returns an S3 class object to allow pretty-printing



# adoptr 0.2.2

* fixed references to SIM paper



# adoptr 0.2.1

* added references to SIM paper in all relevant places



# adoptr 0.2.0

* new feature: composite scores allows more generic expressions than just 
    affine score combinations, cf. composite()
* affine scores (s1 + 2*s2) are no longer supported, use composite instead
* more consistent class system: conditional scores no longer need a specification
    of distributions by default (no need for conditional sampel size e.g.).
    Instead, expected() now requires explicit specification of the data and
    prior distribution to integrate with.
* Vignettes updated
* fixed broked tests due to updated rpact



# adoptr 0.1.1

* extended Description field in DESCRIPTION to full paragraph
* provided examples for all user facing functions
* revision of docs



# adoptr 0.1.0

* initial release
