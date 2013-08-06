\name{krige}
\alias{krige}

\title{
	Kriging function
}
\description{
Variogram analysis and kriging prediction used to analyze feedback and change-in-feedback across space.
}
\usage{
krige(coordinates, statistic, variog.param, grid, krige.param, plots = TRUE)
}

\arguments{
  \item{coordinates}{a 2-column matrix with latitudes and longitudes of observation sites.
}
\item{statistic}{a numeric vector specifying the values, at observation sites, of the statistic to be predicted.
}

\item{variog.param}{a list of arguments used to fit the semivariogram to data, some of these arguments being input of the function \code{variofit} in the \code{geoR} package:
 \describe{
\item{\code{keep.distance}:}{a numeric value between zero and one allowing to specify the maximum distance when fitting the variogram. This maximum distance is the product between the proportion \code{keep.distance} and the maximum distance between observation sites.}

\item{\code{nb.bin}:}{a numeric value specifying the number of bins to be used when fitting the semivariogram.}

\item{\code{cov.pars}:}{a numeric vector with initial values for the exponential covariance parameters, namely the partial sill and the range parameter.}
	
\item{\code{nugget}:}{a numeric value for the nugget parameter. Regarded as a fixed value if \code{fix.nugget = TRUE} or as an initial value for the minimization algorithm if \code{fix.nugget = FALSE}.}
	
\item{\code{fix.nugget}:}{logical, indicating whether the nugget parameter should be regarded as fixed (\code{fix.nugget = TRUE}) or should be estimated (\code{fix.nugget = FALSE}).}

\item{\code{coordinates.scaling}:}{a numeric value by which the matrix \code{coordinates} is multipled when fitting the semivariogram. If \code{coordinates.scaling = 1}, raw coordinates are used.}
    
\item{\code{statistic.scaling}:}{a numeric value by which the statistic to be predicted is multipled when fitting the semivariogram. If \code{statistic.scaling = 1}, the raw statistic is used.}
}
}

  \item{grid}{a list of arguments defining the grid over which the statistic is predicted:
 	\describe{
		\item{\code{x}:}{a numeric vector of longitudes.}
		\item{\code{y}:}{a numeric vector of latitudes.}
		\item{\code{border}:}{a character string that names the region in the world within which the statistic will be predicted. This is a character string in the database \code{world}; see the function \code{map} in the \code{maps} package.}
		\item{\code{proj}:}{a projection definition; see the function \code{project} in the \code{proj4} package.}
		\item{\code{degrees}:}{a logical indicating  if  the lat/long coordinates are assumed to be in degrees (if \code{degrees = TRUE}), or in radians (if \code{degrees = FALSE}).}
	}
	}

  \item{krige.param}{
a list of arguments specifying the type of kriging to be performed and that are used as arguments of the function \code{krige.control} in the \code{geoR} package:
 	\describe{
		\item{\code{type.krige}:}{a character string specifying the type of kriging to be performed. Options are \code{"SK"} and \code{"OK"} corresponding to simple or ordinary kriging. Kriging with external trend and universal kriging can be defined setting \code{type.krige = "OK"} and specifying the trend model using the arguments \code{trend.d}.}

		\item{\code{trend.d}:}{a character string specifying the trend (covariate) values at the data locations.} 
	}
}
  \item{plots}{
a logical indicating if plots characterizing the variogram analysis and the kriging prediction are produced (if \code{plots = TRUE}) or not (if \code{plots = FALSE}), default is TRUE.
}
}
\details{
	If \code{plots = TRUE}, six plots are produced. Plot 1: Estimation of the semivariogram of the statistic using the weighted least squares method (dots: sample semivariogram; curve: theoretical semivariogram -- exponential covariance function with nugget effect). Plot 2: Assessment of the existence of a spatial structure based on a Monte-Carlo test. Plot 3: Assessment of the confidence envelopes of the semivariogram. Plot 4: Boxplots of kriging prediction (left) and kriging standard error (right). Plot 5: Kriging prediction. Plot 6: Kriging standard error.
}
\value{
a list of items characterizing the variogram analysis and the kriging prediction:
\item{input}{the list of arguments in the call of the \code{krige} function except the argument \code{plots}.
}
\item{MAP}{a list allowing to draw the border of the study region that can be made of several polygons. This is the output of the function \code{map} in the \code{maps} package.}
\item{grid}{a 2-column matrix providing the coordinates (in degrees) of the nodes of the prediction locations.}
\item{in.region}{a logical vector indicating, for each grid node whose coordinates are given in  the 2-column matrix \code{grid}, if the node is in the study region or not. The study region is defined by \code{MAP}.}
\item{variofit.wls}{a list providing the result of the weighted-least-square estimation of the semivariogram parameters. This is the output of the function \code{variofit} in the \code{geoR} package.}
\item{krige}{a list providing the result of the spatial prediction for fixed covariance parameters given by \code{variofit.wls}. This is the output of the function \code{krige.conv} in the \code{geoR} package.}
}
\references{
Soubeyrand, S., Morris, C. E. and Bigg, E. K. (2013). Analysis of fragmented time directionality in time series to elucidate feedbacks in climate data. Research Report.
}
\author{
	Samuel Soubeyrand \email{Samuel.Soubeyrand@avignon.inra.fr}, Cindy E. Morris, E. Keith Bigg.
}
\note{
	The \code{krige} function only uses some of the functionalities of the \code{map} function in the \code{maps} package, and the \code{variofit.wls}, \code{krige.conv} and \code{krige.control} functions in the \code{geoR} package. To fully exploit the functionalities of these external functions, the \code{krige} function should be extended.
}


\seealso{\code{rain.feedback.stats}, \code{map} in the \code{maps} package, \code{variofit.wls}  and \code{krige.conv} in the \code{geoR} package
}
\examples{
#### load data of feedback and change-in-feedback indices in 88 sites across Australia
data(rain.feedback.stats)

#### spatial coordinates of the 88 sites
coord=rain.feedback.stats[,3:4]

#### feedback index
stat1=rain.feedback.stats[["Feedback.whole.period"]]

#### variogram analysis and kriging of feedback index
#uncomment:
#par(mfrow=c(2,3), mar=c(5.1,4.1,4.1,4.1))
#kr1=krige(coordinates=coord, statistic=stat1,
#   grid=list(x=seq(110,155,0.25),y=seq(-45,-11,0.25),border="Australia", 
#      proj="+proj=lcc +lat_1=-18 +lat_2=-36 +lat0=-25 +lon_0=140",degrees=TRUE),
#   variog.param=list(keep.distance=2/3,nb.bin=15,cov.pars=c(0.0008,10^6),
#      nugget=0.0001,fix.nugget=FALSE,coordinates.scaling=10^5,statistic.scaling=0.01),
#   krige.param=list(type.krige="ok",trend.d="1st"),plots=TRUE)
}

\keyword{ misc }
