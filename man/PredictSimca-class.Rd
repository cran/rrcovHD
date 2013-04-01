\name{PredictSimca-class}
\Rdversion{1.1}
\docType{class}
\alias{PredictSimca-class}
\alias{show,PredictSimca-method}

\title{Class \code{"PredictSimca"} - prediction of \code{"Simca"} objects}
\description{ The prediction of a "Simca" object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PredictSimca", ...)}
but most often by invoking \code{predict()} on a \code{"Simca"} object. 
They contain values meant for printing by \code{show()}}
\section{Slots}{
  \describe{
    \item{\code{classification}:}{Object of class \code{"factor"} ~~ }
    \item{\code{odsc}:}{A \code{"matrix"} containing the standartized orthogonal distances for each group}
    \item{\code{sdsc}:}{A \code{"matrix"} containing the standartized score distances for each group}
    \item{\code{ct}:}{re-classification table of the training sample}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "PredictSimca")}: Prints the results.. }
	 }
}
\references{
  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\seealso{
  \code{\link{Simca-class}}
}
\examples{
showClass("PredictSimca")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

