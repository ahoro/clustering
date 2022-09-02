\name{InitSoftKMeansModel}
\alias{InitSoftKMeansModel}
\alias{InitSoftKMeansModelFromData}
\title{
Soft K Means Model
}
\description{
Sets up a Soft K Means model for fitting to some data.
}
\usage{
model <- InitSoftKMeansModel(means, beta)
model <- InitSoftKMeansModelFromData(data, k, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{means}{
  The initial means, which will be adjusted as the model iterates, as a d by k matrix,
  where d is the number of data dimensions.
  }
  \item{data}{
  The data that will be fitted, as an n by d matrix, where n is he number of data points.
  }
  \item{k}{
  The number of clusters to fit to the data.
  }
  \item{beta}{
  A stiffness parameter. A lower beta allows more overlap between clusters.
  }
}
\details{
  These functions only initialize the model, they do not perform any fitting. InitSoftKMeansModelFromData
  will choose initial means/clusters uniformly from the range of data supplied.
}
\value{
A SoftKMeansModel object, to be passed to StepModel().
}
\author{
  Adrian Horodeckyj
}