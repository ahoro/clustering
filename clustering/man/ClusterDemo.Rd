\name{ClusterDemo}
\alias{ClusterDemo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering Demonstration
}
\description{
Demonstrates performing clustering with the ArbitraryGMM object. 
}
\usage{
ClusterDemo()
}
%- maybe also 'usage' for other objects documented here.
\details{
The function generates some random data from some random parameters. It then fits 
a gaussian mixture model by EM to it to try to recover those parameters. These are
shown graphically using the rgl package. The generating clusters are shown as
red ellipsoids. The recovered clusters are shown as black wire frame ellipsoids.
The change in means of the fitted gaussians step to step are shown as series of 
arrows.
}

\author{
Adrian Horodeckyj
}
