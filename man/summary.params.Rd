% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/params.R, R/params_lm.R, R/params_mlogit.R,
%   R/params_mlogit_list.R, R/params_surv.R, R/params_surv_list.R
\name{summary.params}
\alias{summary.params}
\alias{summary.params_lm}
\alias{summary.params_mlogit}
\alias{summary.params_mlogit_list}
\alias{summary.params_surv}
\alias{summary.params_surv_list}
\title{Summarize parameter objects}
\usage{
\method{summary}{params_lm}(object, probs = c(0.025, 0.975), ...)

\method{summary}{params_mlogit}(object, probs = c(0.025, 0.975), ...)

\method{summary}{params_mlogit_list}(object, probs = c(0.025, 0.975), ...)

\method{summary}{params_surv}(object, probs = c(0.025, 0.975), ...)

\method{summary}{params_surv_list}(object, probs = c(0.025, 0.975), ...)
}
\arguments{
\item{object}{An object of the appropriate class.}

\item{probs}{A numeric vector of probabilities with values in \verb{[0,1]} used
to compute quantiles. By default, the 2.5th and 97.5th percentiles are
computed.}

\item{...}{Additional arguments affecting the summary. Currently unused.}
}
\value{
A \code{\link[data.table:data.table]{data.table::data.table}} that always contains the following columns:
\describe{
\item{term}{The regression term.}
\item{mean}{The mean value of the regression term.}
\item{sd}{The standard deviation of the values of the regression term.}
}

In addition, the \code{probs} argument determines the quantiles that are computed.
By default, the columns \verb{2.5\%} and \verb{97.5\%} are returned corresponding to the
2.5th and 97.5th percentiles.

Finally, the following columns may also be present:
\describe{
\item{parameter}{The name of the parameter of interest. This is relevant
for any parametric model in which the underlying probability distribution
has multiple parameters. For instance, both \code{\link{params_surv}} and \code{\link{params_surv_list}}
store regression coefficients that are used to model the underlying parameters
of the survival distribution (e.g., shape and scale for a Weibull model). Similarly,
there are two parameters (\code{mean} and \code{sd}) for \code{\link{params_lm}} objects.}
\item{model}{The name of the statistical model. This is used for a
\code{\link{params_surv_list}} object, where each list element represents a separate model.
In a state transition model, each model is a unique health state transition and
in a partitioned survival model, there is a separate model for each curve.}
\item{to}{The health state that is being transitioned to. In \code{\link{params_mlogit}}
and \code{\link{params_mlogit_list}} objects, there are coefficients for each health
state that can be transitioned to.}
\item{from}{The health state that is being transitions from. This is used
for a \code{\link{params_mlogit_list}} objects where a different multinomial
logistic regression is used for each state that can be transitioned from.}
}
}
\description{
Summarize the coefficients of a parameter object by computing the mean,
standard deviation, and quantiles for each model term.
This is a convenient way to check whether a parameter object has been specified
correctly and sampling distributions of the coefficients are as expected.
}
\seealso{
For examples, see the the underlying parameter object functions:
\code{\link[=params_lm]{params_lm()}}, \code{\link[=params_surv]{params_surv()}}, \code{\link[=params_surv_list]{params_surv_list()}}, \code{\link[=params_mlogit]{params_mlogit()}}, and
\code{\link[=params_mlogit_list]{params_mlogit_list()}}.
}
