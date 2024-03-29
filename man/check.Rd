% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hesim_data.R, R/input_mats.R, R/params_lm.R,
%   R/params_mlogit.R, R/params_surv.R, R/params_surv_list.R, R/tparams_mean.R,
%   R/tparams_transprobs.R, R/utils.R
\name{check.id_attributes}
\alias{check.id_attributes}
\alias{check.input_mats}
\alias{check.params_lm}
\alias{check.params_mlogit}
\alias{check.params_surv}
\alias{check.params_surv_list}
\alias{check.tparams_mean}
\alias{check.tparams_transprobs}
\alias{check}
\title{Input validation for class objects}
\usage{
\method{check}{id_attributes}(object, ...)

\method{check}{input_mats}(object, ...)

\method{check}{params_lm}(object, ...)

\method{check}{params_mlogit}(object, ...)

\method{check}{params_surv}(object, ...)

\method{check}{params_surv_list}(object, ...)

\method{check}{tparams_mean}(object, ...)

\method{check}{tparams_transprobs}(object, ...)

check(object, ...)
}
\arguments{
\item{object}{object to check.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
If validation is successful, returns the object in question; otherwise,
informs the user that an error has occurred.
}
\description{
\code{check} is a generic function for validating the inputs of class objects.
}
\keyword{internal}
