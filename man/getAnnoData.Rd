% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determineArrayType.R
\name{getAnnoData}
\alias{getAnnoData}
\title{Get Annotation Object}
\usage{
getAnnoData(type)
}
\arguments{
\item{type}{can be "450k" or "EPIC"}
}
\value{
annotation object, retrieved via getAnnotation()
for the respective array design
}
\description{
Returns the getAnnotation() object for 450k or EPIC data
}
\examples{
getAnnoData("EPIC")[1:3,]
}

