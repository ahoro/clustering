\name{UpgradeModel}
\alias{UpgradeModel}
\title{
Upgrade Model
}
\description{
UpgradeModel can convert a SphericalGMM object to an AxisAlignedGMM object or 
ArbitraryGMM object.
}
\usage{
UpgradeModel(model, modelType)
}
\arguments{
  \item{model}{
  Model object.
  }
  \item{modelType}{
  One of: "AxisAlignedGMM" or "ArbitraryGMM"
  }
}

\value{
  The upgraded object.
}

\author{Adrian Horodeckyj}