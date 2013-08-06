\name{KT.output-class}
\Rdversion{1.1}
\docType{class}
\alias{KT.output-class}
\alias{KT.output}
\alias{KT.output class}
\alias{KT.output-object}
\alias{KT.output object}
\alias{[,KT.output-method}
\alias{[<-,KT.output-method}
\alias{show,KT.output-method}
\alias{summary,KT.output-method}
\alias{names,KT.output-method}
\alias{plot,KT.output-method}

\title{Class \code{"KT.output"}}
\description{
Output of the \link{krige.test} function.
}
\section{Objects from the Class}{
Objects are created by calls of the \link{krige.test} function.
}
\section{Slots}{
  \describe{
    \item{\code{krige.output}:}{Object of class \code{"list"}, output of the \link{krige} function.}
    \item{\code{subregion}:}{Object of class \code{"list"}, two vectors \code{x} and \code{y} containing the latitudes and the longitudes, respectively, of the vertices of a polygon. The polygon defines a subregion where one supposes variation in the predicted index.}
    \item{\code{averageKrigingPrediction.rand}:}{Object of class \code{"numeric"} specifying the averages of the kriging predictions in \code{subregion} obtained with randomized data (here, a randomization is a random translation on a torus). }
    \item{\code{averageKrigingPrediction.obs}:}{Object of class \code{"numeric"} specifying the average of the kriging prediction in \code{subregion} obtained with observed data. }
    \item{\code{alternative}:}{Object of class \code{"character"}, "greater" or "less".}
    \item{\code{p.value}:}{Object of class \code{"numeric"}, p-value of the test.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "KT.output")}}
    \item{[<-}{\code{signature(x = "KT.output")}}
    \item{show}{\code{signature(object = "KT.output")}}
    \item{summary}{\code{signature(object = "KT.output")}}
    \item{names}{\code{signature(object = "KT.output")}}
    \item{plot}{\code{signature(object = "KT.output")}}
	 }
}
\references{
Soubeyrand, S., Morris, C. E. and Bigg, E. K. (2013). Analysis of fragmented time directionality in time series to elucidate feedbacks in climate data. Research Report.
}
\author{
	Samuel Soubeyrand \email{Samuel.Soubeyrand@avignon.inra.fr}, Cindy E. Morris, E. Keith Bigg.
}
\examples{
showClass("KT.output")

#### load data of feedback and change-in-feedback indices in 88 sites across Australia
data(rain.feedback.stats)

#### spatial coordinates of the 88 sites
coord=rain.feedback.stats[,3:4]

#### map of feedback index computed from the whole data series
stat1=rain.feedback.stats[["Feedback.whole.period"]]

#### variogram analysis and kriging of feedback index
#uncomment:
#par(mfrow=c(2,3), mar=c(5.1,4.1,4.1,2.1))
#kr1=krige(coordinates=coord,statistic=stat1,
#   grid=list(x=seq(110,155,0.25),y=seq(-45,-11,0.25),border="Australia",
#      proj="+proj=lcc +lat_1=-18 +lat_2=-36 +lat0=-25 +lon_0=140",degrees=TRUE),
#   variog.param=list(keep.distance=2/3,nb.bin=15,cov.pars=c(0.0008,10^6),
#      nugget=0.0001,fix.nugget=FALSE,coordinates.scaling=10^5,statistic.scaling=0.01),
#   krige.param=list(type.krige="ok",trend.d="1st"),plots=TRUE)

#### test spatial variation in feedback index
## computer intensive stage (un-comment the following command lines)
#uncomment:
#kt1=krige.test(krige.output=kr1,subregion=list(x=c(138,152,152,138),y=-c(40,40,33,33)),
#   alternative="greater", nb.rand=2000)
	
## plot test output	
#uncomment:
#par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,2.1))	
#plot(kt1,digits=list(predict=3,pvalue=3),breaks=12)
}
\keyword{classes}