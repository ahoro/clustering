\name{StepModel}
\alias{StepModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting Step for Models
}
\description{
StepModel does one iteration of the relevant clustering algorithm (K means or EM). It is 
up to the caller to decide when to stop. This can be done, for example, either by running
for a predetermined number of steps or by checking whether the model parameters are no 
longer significantly changing, for some definition of "significantly".
}
\usage{
model <- StepModel(model, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
An object returned by one of the Init functions.
}
  \item{data}{
Data for the model to be fitted to, as a n by d matrix, where n is the number of data points
and d is the number of dimensions of the data.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The returned object will be a model object of the same type with updated model parameters.
}
\references{
  D. J. C. MacKay, "Information Theory, Inference, and Learning Algorithms", https://www.inference.org.uk/mackay/itila/ Chapters 20 & 22
  
  Bill Press, "Opinionated Lessons in Statistics", http://opinionatedlessons.org/ Segment 20
  
}
\author{
  Adrian Horodeckyj
}