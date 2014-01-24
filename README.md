GAMMKS
========================================================

R package source for developing Generalized Addititive Multiple Models with Kalman Smoothing (GAMMKS) for load estimtion.  Load estimates are based on irregular or infrequent water quality measurements and daily streamflow data.  GAMMKS estimates a generalized additive model (GAM) and a linear model (LM) on the basis of explanatory variables specified by the user. In general, the GAM model provides greater flexibility for modeling the relation between the explanatory and response variables when data are available to support estimation and predicition.  The LM provides greater robustness for estimation during periods when their are gaps in water-quality data collection or the water-quality data is sparse.  Both models require daily streamflow information to faciliate estimation.  

Once the GAM and LM are estimated, they are used to estimate the magnitudes and uncertainties of daily concentrations across the target period of record using explanatory variables that commonly describe depenence on: 

1. Daily streamflow

2. Seasonal components based on date information

3. Trend components based on date information

4. Optionally using flow anomaly data [**waterData** package in R.](http://cran.r-project.org/)

Daily estimtes of concentrations computed by the two models are combined inversely proportional to their respective variances.  The integrated estimate is referred to as the GAMM estimates.

A Kalman Smoother is then used to integrate information from direct water quality measurements with the GAMM estimate.


