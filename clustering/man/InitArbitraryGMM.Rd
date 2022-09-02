\name{InitArbitraryGMM}
\alias{InitArbitraryGMM}
\alias{InitArbitraryGMMFromData}
\title{
Arbitrary GMM
}
\description{
Sets up a gaussian mixture model to be fitted to some data with EM. Each cluster
has a covariance matrix.
}
\usage{
model <- InitArbitraryGMM(means, sigmas, priors)
model <- InitAxisAlignedGMMFromData(data, k)
}
\arguments{
  \item{means}{
  The initial means, which will be adjusted as the model iterates, as a d by k matrix,
  where d is the number of data dimensions.
  }
  \item{sigmas}{
  The initial covariance matrices, which will be adjusted as the model iterates, a list
  of length k of d by d matrices.
  }
  \item{priors}{
  The prior probability of a data point belonging to a particular cluster. In other words,
  the proportion of points that belong to a given cluster. A vector of length k.
  }
  \item{data}{
  The data that will be fitted, as an n by d matrix, where n is he number of data points.
  }
  \item{k}{
  The number of clusters to fit to the data.
  }
}
\details{
  These functions only initialize the model, they do not perform any fitting. InitArbitraryGMMFromData
  will choose initial means/clusters uniformly from the range of data supplied, standard deviations of 1, 
  aligned to the axes, and uniform priors.
}
\value{
  An ArbitraryGMM object, to be passed to StepModel()
}
\author{
  Adrian Horodeckyj
}
