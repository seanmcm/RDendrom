# RDendrom
An R package for analysis, diagnostics, and presentation of intra-annual tree growth time series.

\code{RDendrom} is a suite of functions designed to handle time series data from 'manual' dendrometer band measurements conducted intensively through the growing season. The data enter in a set 'long' form with such as in the following header file.

[instert header of a data.frame]

The general workflow for these data includes processing, analyzing, vetting, and presenting individual annual tree growth. There are also functions to aggregate results, which are more generally in 'beta' (they have been used for particular analyses, but not formalized or generalized for use in the package yet, and so need a bit of attention).

## Processing

The primary approach to 'processing' data is through fitting a three-part general additive model (GAM) to the data. This generally follows McMahon and Parker 2014, Ecology and Evolution. The core of the procedure is the fitting of a 5-parameter logistic function to the growth time series. The two other parts of the GAM are the starting diameter and ending diameter which essentially estimates the leaf-off dormant diamter at breast height (DBH) from the winter before and after the growing season. The LG5 function relates dbh, or the diameter measurements at a given day of the year.  for a designated series of days.

![](LG5.png)

The LG5 comes from a family of functions, known collectively
as Richards models, that can produce good candidate parameterizations for intra-
annual tree growth data \citep{Oswald:2012ux}. We found the effort to fit so many parameters worthwhile
because the parameter number is important for precise fits, necessary to
discover small fluctuations in growth during the year.  Four of the parameters
have straightforward interpretation: two parameters define the lower and upper
asymptotes (\emph{L} and \emph{K} respectively), \emph{doy_{ip}} marks the day of the year
when the inflection of the curve is predicted to occur, and the rate parameter
r describes the slope of the curve at the inflection point. The final
parameter, \theta allows asymmetrical fits by changing the approach of to
the upper asymptote. This is a critical parameter as when \theta = 1, the
curve is symmetrical but growth forms are rarely symmetrical and so its
inclusion merits the increased efforts in fitting a 5-parameter model.  

This fit is estimated using two passes thorugh the \emph{optim()} function.  The resulting curve (which is fit without constraining \emph{L} or \emph{K} to starting and stopping diamters) is then combined with estimates of these diameters of growth initiation and cessation (a and b respectively). The horizontal lines at a and b are then combined with the LG5 curve to produce the GAM that best fits the annual time series data. These three pieces of the GAM are combined in order to allow optimal fits to the growth curve while capturing the sometimes abrupt beginning of bole growth and smoother asymptote of ending of growth. 

## Analyzing

After the initial fits are conducted, additional metrics can be estimated. The primary function for this is the 

## Vetting

The two vetting functions, \emph{id.outliers()} and \emph{find.outliers()}, allow either manual or automatic designation of problem measurements respectively. \emph{id.outliers()} allows interactive (point and click) identification of outlies that are flagged, re-incorporated into the data, and then ignored during re-runs of parameter estimation. \emph{find.outliers()} uses residuals from the fit curves to flag values outside of a designated numer of standard deviations from zero (default = 3). 

It is worth noting that within the \emphs{get.optimized.dendro()} function there are several filters that also serve to vet data. Insufficient sample size (default = 9 values) results in skipping the fitting function. That function also fits a linear model to the data and if the slope is negative (representing a dead or shrinking stem), then fitting is skipped. The fits also produce R^2 values for every fit that is implemented, so that further vetting and subsetting of data for specific analyses can take into account how well the curves fit the data.

## Presenting

There are several plotting functions that present the data directly, with fits, with hulls, over years, and aggregated. Help files explain these different plotting possibilities.
