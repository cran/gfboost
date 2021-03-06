% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RejStep.R
\name{RejStep}
\alias{RejStep}
\title{CMB validation step}
\usage{
RejStep(
  D,
  nsing,
  Bsing = 1,
  ind,
  sing = FALSE,
  singfam = Gaussian(),
  evalfam = Gaussian(),
  M = 10,
  m_iter = 100,
  kap = 0.1,
  LS = FALSE,
  best = 1
)
}
\arguments{
\item{D}{Data matrix. Has to be an \eqn{n \times (p+1)-}dimensional data frame in the format \eqn{(X,Y)}. The \eqn{X-}part must not
contain an intercept column containing only ones since this column will be added automatically.}

\item{nsing}{Number of observations (rows) used for the SingBoost submodels.}

\item{Bsing}{Number of subsamples based on which the SingBoost models are validated. Default is 1. Not to confuse with parameter \code{B} for the Stability Selection.}

\item{ind}{Vector with indices for dividing the data set into training and validation data.}

\item{sing}{If \code{sing=FALSE} and the \code{singfam} family is a standard Boosting family that is contained in the package
\code{mboost}, the CMB aggregation procedure is executed for the corresponding standard Boosting models.}

\item{singfam}{A SingBoost family. The SingBoost models are trained based on the corresponding loss function. Default is \code{Gaussian()} (squared loss).}

\item{evalfam}{A SingBoost family. The SingBoost models are validated according to the corresponding loss function. Default is \code{Gaussian()} (squared loss).}

\item{M}{An integer between 2 and \code{m_iter}. Indicates that in every \eqn{M-}th iteration, a singular iteration will be
performed. Default is 10.}

\item{m_iter}{Number of SingBoost iterations. Default is 100.}

\item{kap}{Learning rate (step size). Must be a real number in \eqn{]0,1]}. Default is 0.1 It is recommended to use
a value smaller than 0.5.}

\item{LS}{If a \code{singfamily} object that is already provided by \code{mboost} is used, the respective Boosting algorithm
will be performed in the singular iterations if \code{Ls} is set to \code{TRUE}. Default is \code{FALSE}.}

\item{best}{Needed in the case of localized ranking. The parameter \code{K} of the localized ranking loss will be
computed by \eqn{best \cdot n} (rounded to the next larger integer). Warning: If a parameter \code{K} is inserted into the
\code{LocRank} family, it will be ignored when executing SingBoost.}
}
\value{
\item{loss}{Vector of validation losses.}
\item{occ}{Selection frequencies for each Boosting model.}
}
\description{
{Validation step to combine different SingBoost models.}
}
\details{
{Divides the data set into a training and a validation set. The SingBoost models are computed on the training
set and evaluated on the validation set based on the loss function corresponding to the selected Boosting family.}
}
\references{
{Werner, T., Gradient-Free Gradient Boosting, PhD Thesis, Carl von Ossietzky University Oldenburg, 2020}
}
