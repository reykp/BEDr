\name{duos}
\alias{duos}
\docType{package}
\title{
This function perform Bayesian density estimation.
}
\description{
The Bayesian density estimate is attained using the DUOS (Density on Uniform Order Statistics) method.
}
\usage{
duos(y, MH_N=20000, k=12, alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A univariate numeric vector for estimating the density
}
\item{MH_N}{
The number of iterations to perform in the Gibb's sampler
}
\item{k}{
The number of cut-points to use in the algorithm
}
\item{alpha}{
The prior parameter to use in the Dirichlet prior (it is constant)
}
}
\value{
\item{C}{
  A matrix containing the posterior draws for the cut-point parameters
}
\item{P}{
  A matrix containing the poserior draws for the bin probabilities
}
}
