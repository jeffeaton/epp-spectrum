source("spectrum.R")
source("../analysis-functions.R")

theta <- c(-27.2776051, 0.2286067, 1.5936094, 2.0196487, -0.4601538, -1.6417664, 1.4846658, 0.1725820, 0.3623600)

t(matrix(fnBSpline(theta[3:9]), 10))
