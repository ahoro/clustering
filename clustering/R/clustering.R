# Copyright 2022 Adrian Horodeckyj

# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

StepModel <- function(model, data) {
  if(model$d != dim(data)[2]) {
    stop("Number of dimensions in the model does not match the data.")
  }
  UseMethod("StepModel", model)
}

InitialModel <- function(means) {
  if(!is.matrix(means)) {
    stop("Input must be a d by k matrix.")
  }
  model <- list()
  model$d = dim(means)[1]
  model$k = dim(means)[2]
  model$means = means
  return(model)
}

InitialModelFromData <- function(data, k) {
  dataMins <- apply(X = data, FUN = min, MARGIN = 2)
  dataMaxes <- apply(X = data, FUN = max, MARGIN = 2)
  model <- InitialModel(replicate(n = k, runif(n = model$d, min = dataMins, max = dataMaxes))) # d by k
  return(model)
} 

# Create initial model by specifying the starting means directly as a d by k matrix
InitHardKMeansModel <- function(means) {
  model <- InitialModel(means)
  class(model) <- "HardKMeansModel"
  return(model)
}

# Create initial model from the dimensionality of the data and given number of means 
InitHardKMeansModelFromData <- function(data, k) {
  model <- InitialModelFromData(data, k)
  class(model) <- "HardKMeansModel"
  return(model)
}

# "Hard" K means algorithm, Algorithm 20.2 in MacKay.
StepModel.HardKMeansModel <- function(model, data) {
  ind <- apply(X = data,  MARGIN = 1, 
               FUN = function(x) {
                 return(which.min(colSums((model$means - x)^2)))
               })
  newMeans <- sapply(X = split(data, factor(ind, 1:model$k)), 
                     FUN = function(cluster) {
                       return(rowMeans(matrix(data = cluster, nrow = model$d)))
                     })
  nonEmptyClusters <- !is.nan(newMeans[1,]) # rowMeans() creates NaNs where a mean owns no points, which keep old values
  model$means[, nonEmptyClusters] <- newMeans[, nonEmptyClusters] 
  return(model)
}

# Create initial model by specifying the starting means directly as a d by k matrix
InitSoftKMeansModel <- function(means, beta) {
  model <- InitialModel(means)
  model$beta <- beta
  class(model) <- "SoftKMeansModel"
  return(model)
}

# Create initial model from the dimensionality of the data and given number of means 
InitSoftKMeansModelFromData <- function(data, k, beta) {
  model <- InitialModelFromData(data, k)
  model$beta <- beta
  class(model) <- "SoftKMeansModel"
  return(model)
}

# "Soft" K means algorithm, Algorithm 20.7 in MacKay
StepModel.SoftKMeansModel <- function(model, data) {
  resp <- apply(X = model$means, 
                FUN = function(mu) {
                  return(exp(- model$beta * sqrt(colSums(mu - t(data))^2)))
                }, MARGIN = 2)
  resp <- resp / rowSums(resp)
  
  model$means <- t(data) %*% resp / colSums(resp)
    
  return(model)
}

InitGMM <- function(means, priors) {
  if(!is.list(means)) {
    stop("means must be a list of length k of vectors.")
  }
  if(!is.vector(priors) || length(priors) != length(means)) {
    stop("priors must be a vector of length k")
  }
  model <- list()
  model$k <- length(means)
  model$d <- length(means[[1]])
  model$means <- means
  model$priors <- priors
  return(model)
}

InitGMMFromData <- function(data, k) {
  dataMins <- apply(X = data, FUN = min, MARGIN = 2)
  dataMaxes <- apply(X = data, FUN = max, MARGIN = 2)
  model <- InitGMM(means = replicate(n = k, runif(n = dim(data)[2], min = dataMins, max = dataMaxes), simplify = FALSE), # d by k
                   priors = rep(1/k, k))
  return(model)
}

InitSphericalGMM <- function(means, stdDevs, priors) {
  model <- InitGMM(means, priors)
  if(!is.vector(stdDevs) || length(stdDevs) != model$k) {
    stop("stdDevs must be a vector of length k")
  }
  model$sigmaSq <- stdDevs
  class(model) <- "SphericalGMM"
  return(model)
}

InitSphericalGMMFromData <- function(data, k) {
  model <- InitGMMFromData(data, k)
  model$sigmaSq <- replicate(k, 1) # arbitrarily start with std. dev. of 1
  class(model) <- "SphericalGMM"
  return(model)
}

# EM on a GMM with all standard deviations all equal, Algorithm 22.2 in MacKay
# This differs from the book because it uses log responsibilities to avoid 
# numerical underflow.
StepModel.SphericalGMM <- function(model, data) {
  # E step
  dists <- sapply(X = model$means, FUN = function(mu) { return(rowSums(t(mu - t(data))^2)) })
  
  logResp <- log(model$priors / (sqrt(2 * pi * model$sigmaSq))^model$d) - t(t(dists) / model$sigmaSq)
  logResp <- logResp - apply(X = logResp, FUN = LogSumExp, MARGIN = 1)
  
  # M step
  resp <- exp(logResp)
  
  respSumOvern <- colSums(resp)
  
  means <- Rows2List(t(resp) %*% data / respSumOvern)
  
  model$sigmaSq <- mapply(FUN = function(mu, resps) { 
    return(sum(t(mu - t(data))^2 * resps))
  }, model$means, Columns2List(resp)) / respSumOvern
  
  model$priors <- respSumOvern / nrow(data)
  
  return(model)
}

InitAxisAlignedGMM <- function(means, stdDevs, priors) {
  model <- InitGMM(means, priors)
  if(!is.list(stdDevs) 
     || ! Reduce(f = function(b, v) { b && is.vector(v) && length(v) == model$d },
                 x = stdDevs, init = length(sigmas) == model$k)) {
    stop("stdDevs must be a list of k vectors of length d")
  }
  model$sigmas <- sigmas
  class(model) <- "AxisAlignedGMM"
  return(model)
}

InitAxisAlignedGMMFromData <- function(data, k) {
  model <- InitGMMFromData(data, k)
  model$sigmaSq <- replicate(k, rep(1, d), simplify = FALSE) # start with std. dev. of 1 in every axis
  class(model) <- "AxisAlignedGMM"
  return(model)
}

# EM on a GMM with axis-aligned gaussians, with separate standard deviations per dimension. 
# Algorithm 22.4 in MacKay, but calculating with log responsibilities
StepModel.AxisAlignedGMM <- function(model, data) {
  # E step
  dists <- mapply(FUN = function(mu, sigma) {
                          return(- rowSums(t(mu - t(data))^2 / 2 * sigma)) 
                        }, model$means, model$sigmaSq)
  
  coeff <- log(priors / (sqrt(2 * pi * sapply(X = model$sigmaSq, FUN = prod)))) # could change to sum
  
  logResp <- dists + coeff 
  
  logRespSumOverk <- apply(X = logResp, MARGIN = 1, FUN = LogSumExp)
  
  logResp <- logResp - logRespSumOverk
  # M step
  resp <- exp(logResp) # n by k
  
  respSumOvern <- colSums(resp)
  
  model$means <- Rows2List((t(resp) %*% data) / respSumOvern)
  
  model$sigmaSq <- mapply(FUN = function(mu, r, rson) {
                                   return(colSums(r * t(t(data) - mu)^2 / rson))
                                }, model$means, Columns2List(resp), respSumOvern, SIMPLIFY = FALSE)
  
  model$priors <- respSumOvern / nrow(data)
  
  return(model)
}

# log sum of logarithms, helps avoid underflow
LogSumExp <- function(z) {
  mz <- max(z)
  return(mz + log(sum(exp(z - mz))))
}

Rows2List <- function(X) {
  return(apply(X = X, MARGIN = 1, FUN = identity, simplify = FALSE))
}

Columns2List <- function(X) {
  return(apply(X = X, MARGIN = 2, FUN = identity, simplify = FALSE))  
}

InitArbitraryGMM <- function(means, sigmas, priors) {
  model <- InitGMM(means, priors)
  if(!is.list(sigmas) 
     || ! Reduce(f = function(b, m) { b && is.matrix(m) && all(dim(m) == c(model$d, model$d)) },
                 x = sigmas, init = length(sigmas) == model$k)) {
    stop("sigmas must be a list of k d by d matrices")
  }
  model$sigmas <- sigmas
  class(model) <- "ArbitraryGMM"
  return(model)
}

InitArbitraryGMMFromData <- function(data, k) {
  model <- InitGMMFromData(data, k)
  model$sigmas <- replicate(k, diag(model$d), simplify = FALSE) # start with identity matrices
  class(model) <- "ArbitraryGMM"
  return(model)
}

# EM on a GMM with arbitrary gaussians, that is, non-axis-aligned. Exercise 22.7 in MacKay.
# Based on equations from Lecture 28 from Press.
StepModel.ArbitraryGMM <- function(model, data) {
  # E step
  logLiklihoods <- mapply(FUN = function(mu, sigma) {
                                    # Cholesky decomposition is O(n^2) as is backsolve() whereas
                                    # matrix inversion is O(n^3). This is more numerically stable as well.
                                    decSigma <- chol(sigma) # pivot = TRUE?
                                    logDetSigma <- log(det(sigma))
                                    return(apply(X = data, MARGIN = 1, 
                                                 FUN = function(x) {
                                                         dist <- x - mu
                                                         # backsolve() should be numerically more stable than multiplying by an inverted matrix
                                                         return(-0.5 * dist %*% backsolve(decSigma, backsolve(decSigma, dist, transpose = TRUE)) - (model$d / 2) * log(2 * pi) - 0.5 * logDetSigma)
                                                        }))
                                  }, model$means, model$sigmas) # results in n by k matrix
  logResp <- t(t(logLiklihoods) + log(model$priors)) # multiply by prior, in log space, row-wise
  logRespSumOverk <- apply(X = logResp, MARGIN = 1, FUN = LogSumExp) # should be length n
  logResp <- logResp - logRespSumOverk # divide by evidence / normalize, in log space, column-wise
  # M step
  resp <- exp(logResp) # again, n by k
  respSumOvern <- colSums(resp) # length k
  model$means <- Rows2List(X = t(t(data) %*% resp) / respSumOvern)
  # This could have been done by making a 3 dimensional array and summing along one dimension, however
  # this way is clearer because it keeps the matrices as matrices, not array slices. It also prevents
  # having to allocate all the working memory at once. 
  model$sigmas <- mapply(FUN = function(mu, respVec, respSumScalar) {
                                mats <- mapply(FUN = function(x, respScalar) { # a list of matrices
                                                       return(respScalar * ((x - mu) %*% t(x - mu)))
                                                      }, Rows2List(X = data), respVec, SIMPLIFY = FALSE)
                                return(Reduce(f = `+`, x = mats) / respSumScalar)
                              }, model$means, Columns2List(X = resp), respSumOvern, SIMPLIFY = FALSE)
  
  model$priors <- respSumOvern / nrow(data)
  
  return(model)
}

# Allows a GMM to change to a more sophisticated kind of GMM
UpgradeModel <- function(model, modelType) {
  if(class(model) == "SphericalGMM") {
    model <- InitAxisAlignedGMM(means = model$means, 
                                stdDevs = lapply(X = model$sigmaSq, FUN = rep, model$k),
                                priors = model$priors)
  }
  if(class(model) == "AxisAlignedGMM") {
    if(modelType == "AxisAlignedGMM") {
      return(model)
    }
    model <- InitArbitraryGMM(means = model$means, 
                              sigmas = lapply(X = model$sigmaSq, FUN = diag), 
                              priors = model$priors)
  }
  if(class(model) == "ArbitraryGMM" && modelType == "ArbitraryGMM") {
    return(model)
  } 
  stop("Cannot upgrade model to given type.")
}