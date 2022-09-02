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

library(rgl)

# Produce a matrix that represents a rotation by the given angle, in d dimensions, 
# in the plane of dimensions d1, d2
# e.g. A rotation about the x-axis (1) is a rotation in the yz-plane (2, 3)
RotationMatrix <- function(angle, d, d1, d2) {
  if(d1 < 1 || d2 < 1 || d1 == d2 || d1 > d || d2 > d) {
    stop("d1, d2 must be distinct and less than or equal to number of dimensions.")
  }
  m <- diag(x = 1, nrow = d, ncol = d)
  m[d1, d1] <- cos(angle)
  m[d1, d2] <- -sin(angle)
  m[d2, d1] <- sin(angle)
  m[d2, d2] <- cos(angle)
  return(m)
}

ClusterDemo <- function() {
  data <- NULL
  clusters <- sample(2:8, 1) # The number of clusters in the data
  open3d()
  for(i in 1:clusters) {
    # Randomize the parameters of the random data
    n <- 10 * rpois(n = 1, lambda = 20)
    m <- rnorm(n = 3, sd = 2)
    s <- rgamma(n = 3, shape = 1.5, scale = 0.5)
    samples <- matrix(data = rnorm(n = n * 3, mean = 0, sd = s), nrow = n, ncol = 3, byrow = TRUE)
    # Rotate in each dimension by a random angle
    xRot <- RotationMatrix(runif(n = 1, min = 0, max = 2 * pi), 3, 2, 3)
    yRot <- RotationMatrix(runif(n = 1, min = 0, max = 2 * pi), 3, 1, 3)
    zRot <- RotationMatrix(runif(n = 1, min = 0, max = 2 * pi), 3, 1, 2)
    r <- xRot %*% yRot %*% zRot
    samples <- samples %*% r
    samples <- t(t(samples) + m) # shift to chosen means
    covariance <- t(r) %*% diag(s) %*% r # this is how you transform a covariance matrix, note it's not just a rotation
    plot3d(ellipse3d(x = covariance, centre = m, level = 0.68), col = "red", alpha = 0.1, type = "shade", add = TRUE)
    data <- rbind(data, samples)
  }
  plot3d(data, add = TRUE)
  # Assume we know the number of clusters, but in general, this is not the case
  model <- InitArbitraryGMMFromData(data, k = clusters) 
  for(i in 1:10) { # This tends to converge pretty quickly, so just take 10 steps
    newModel <- StepModel(model, data)
    for(j in 1:model$k) {
      arrow3d(model$means[[j]], newModel$means[[j]], type = "lines")
    }
    model <- newModel
  }
  for(i in 1:model$k) {
    plot3d(ellipse3d(x = model$sigmas[[i]], centre = model$means[[i]], level = 0.68), alpha = 0.1, type = "wire", add = TRUE)
  }
}
