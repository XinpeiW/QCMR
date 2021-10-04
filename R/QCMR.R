#' @title Quadratic Mendelian randomization analysis
#
#' @description This package is used to identify and estimate the quadratic causality of exposure on outcome.
#'
#' @param data
#'
#' @return res
#'
#' @examples
#'
#' @export


QCMR <- function(data) {
  data <- na.omit(data)
  Z <- data$Z; X <- data$X; Y <- data$Y

  # coefficients estimation
  obj1 <- summary(lm(X~Z))$coefficients
  obj2 <- summary(lm(Y~Z+I(Z^2)))$coefficients
  ahat <- obj1[2,1]
  beta2hat <- obj2[3,1]/ahat^2
  beta1hat <- (obj2[2,1]-2*ahat*beta2hat*obj1[1,1])/ahat

  # se calculation
  d <- data.frame(cbind(Z,X,Y))
  library(bootstrap)
  f <- function(raw,d){
    objx <- summary(lm(d[raw,2] ~ d[raw,1]))$coefficients
    objy <- summary(lm(d[raw,3] ~ d[raw,1] + I(d[raw,1]^2)))$coefficients
    ah <- objx[2,1]
    beta2h <- objy[3,1]/ah^2
    beta1h <- (objy[2,1]-2*ah*beta2h*objx[1,1])/ah
    return(cbind(beta1h,beta2h))
  }
  result <- bootstrap(1:nrow(d),nboot=100,theta=f,d)
  se1 <- sd(result$thetastar[1,])
  se2 <- sd(result$thetastar[2,])
  p1 <- (1-pnorm(abs(beta1hat[i]/se1[i])))*2
  p2 <- (1-pnorm(abs(beta2hat[i]/se2[i])))*2

  res <- cbind(beta1hat,se1,p1,beta2hat,se2,p2)
  return(res)

}
