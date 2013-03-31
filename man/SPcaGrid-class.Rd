\name{SPcaGrid-class}
\Rdversion{1.1}
\docType{class}
\alias{SPcaGrid-class}
\alias{getQuan,SPcaGrid-method}

\title{Class \code{SPcaGrid} - Sparse Robust PCA using PP - GRID search Algorithm}
\description{
    Holds the results of an approximation of the PP-estimators for sparse and robust PCA using the grid search algorithm in the plane.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SPcaGrid", ...)} but the
    usual way of creating \code{SPcaGrid} objects is a call to the function
    \code{SPcaGrid()} which serves as a constructor.
}



\section{Slots}{
  \describe{
   
    \item{\code{call}, \code{center}, \code{scale}, \code{loadings},
        \code{eigenvalues}, \code{scores}, \code{k},
        \code{sd}, \code{od}, \code{cutoff.sd}, \code{cutoff.od},
        \code{flag}, \code{n.obs}:}{
        from the \code{"\link[rrcov]{Pca-class}"} class.}
    
  }
}
\section{Extends}{
Class \code{"\linkS4class{PcaGrid}"}, directly.
Class \code{"\linkS4class{PcaRobust}"}, by class "PcaGrid", distance 2.
Class \code{"\linkS4class{Pca}"}, by class "PcaGrid", distance 3.
}
\section{Methods}{
  \describe{
    \item{getQuan}{\code{signature(obj = "SPcaGrid")}: ... }
	 }
}
\references{
 Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{SPcaGrid}}, \code{\link{PcaGrid-class}}, \code{\link{PcaRobust-class}}, \code{\link{Pca-class}}, \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}
}
\examples{
showClass("SPcaGrid")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

