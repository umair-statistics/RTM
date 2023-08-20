\name{LeadExposedChildren}
\alias{LeadExposedChildren}
\docType{data}
\title{Blood Lead Levels from the
Treatment of Lead Exposed Children (TLC) Trial.}
\description{The Treatment of Lead-Exposed Children (TLC) trial was a placebo-controlled, randomized study of succimer (a chelating agent) in children with blood lead levels of 20-44 micrograms/dL. These data consist of four repeated measurements of blood lead levels obtained at baseline (or week 0), week 1, week 4, and week 6 on 100 children who were randomly assigned to chelation treatment with succimer or placebo.}
\usage{data("LeadExposedChildren")}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{\code{ID}}{a numeric vector represent identity of the patients (1--100)}
    \item{\code{Group}}{a character represent the treatment group (P=placebo, A=treatment)}
    \item{\code{Week0}}{a baseline measuremnts of Blood lead level}
    \item{\code{Week1}}{a repeated followup measuremnts of Blood lead level after first week}
    \item{\code{Week4}}{a repeated followup measuremnts of Blood lead level after fourth week}
    \item{\code{Week6}}{a repeated followup measuremnts of Blood lead level after sixth week}
  }
}
\source{
Data courtesy of Dr. George G. Rhoads (Chair, TLC Steering Committee).}
\references{
Treatment of Lead-exposed Children (TLC) Trial Group. (2000).
Safety and Efficacy of Succimer in Toddlers with Blood Lead Levels of
20-44 µg/dL. \emph{Pediatric Research}, 48, 593-599.}
\examples{
data(LeadExposedChildren)
str(LeadExposedChildren)
}
\keyword{datasets}
