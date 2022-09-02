\name{InitAxisAlignedGMM}
\alias{InitAxisAlignedGMM}
\alias{InitAxisAlighnedGMMFromData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Axis-Aligned GMM
}
\description{
Sets up an axis-aligned gaussian mixture model to be fitted to some data with EM. The standard
deviation for each cluster is fitted per data dimension.
}
\usage{
model <- InitAxisAlignedGMM(means, stdDevs, priors)
model <- InitAxisAlignedGMMFromData(data, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{means}{
  The initial means, which will be adjusted as the model iterates, as a d by k matrix,
  where d is the number of data dimensions.
  }
  \item{stdDevs}{
  The initial standard deviations, which will be adjusted as the model iterates, as a list
  of length k of vectors of length d.
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
  These functions only initialize the model, they do not perform any fitting. InitAxisAlignedGMMFromData
  will choose initial means/clusters uniformly from the range of data supplied, standard deviations of 1,
  and uniform priors.
}
\value{
  An AxisAlignedGMM object, to be passed to StepModel()
}
\author{
  Adrian Horodeckyj
}