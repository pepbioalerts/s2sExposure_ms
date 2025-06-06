# TOPSIS vector ----
TOPSISv_2 <- function (decision, weights, cb) 
{
  if (!is.matrix(decision)) 
    stop("'decision' must be a matrix with the values of the alternatives")
  if (missing(weights)) 
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if (!all.equal(sum(weights),1))
    stop("The sum of 'weights' is not equal to 1")
  if (!is.character(cb)) 
    stop("'cb' must be a character vector with the type of the criteria")
  if (!all(cb == "max" | cb == "min")) 
    stop("'cb' should contain only 'max' or 'min'")
  if (length(weights) != ncol(decision)) 
    stop("length of 'weights' does not match the number of the criteria")
  if (length(cb) != ncol(decision)) 
    stop("length of 'cb' does not match the number of the criteria")
  d = sqrt(colSums(decision^2))
  NW <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    NW[, j] <- (decision[, j]/d[j]) * weights[j]
  }
  posI <- as.integer(cb == "max") * apply(NW, 2, max) + as.integer(cb == 
                                                                     "min") * apply(NW, 2, min)
  negI <- as.integer(cb == "min") * apply(NW, 2, max) + as.integer(cb == 
                                                                     "max") * apply(NW, 2, min)
  distance = function(x, y) {
    sqrt(sum((x - y)^2))
  }
  posDis <- apply(NW, 1, distance, posI)
  negDis <- apply(NW, 1, distance, negI)
  R <- negDis/(negDis + posDis)
  return(data.frame(Alternatives = 1:nrow(decision), R = R, 
                    Ranking = rank(-R, ties.method = "first")))
}
# TOPSIS linear ----
TOPSISl_2 <- function (decision, weights, cb) 
{
  if (!is.matrix(decision)) 
    stop("'decision' must be a matrix with the values of the alternatives")
  if (missing(weights)) 
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if (!all.equal(sum(weights),1))
    stop("The sum of 'weights' is not equal to 1")
  if (!is.character(cb)) 
    stop("'cb' must be a character vector with the type of the criteria")
  if (!all(cb == "max" | cb == "min")) 
    stop("'cb' should contain only 'max' or 'min'")
  if (length(weights) != ncol(decision)) 
    stop("length of 'weights' does not match the number of the criteria")
  if (length(cb) != ncol(decision)) 
    stop("length of 'cb' does not match the number of the criteria")
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  Norm <- as.integer(cb == "max") * apply(decision, 2, max) + 
    as.integer(cb == "min") * apply(decision, 2, min)
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    if (cb[j] == "max") {
      N[, j] <- decision[, j]/Norm[j]
    }
    else {
      N[, j] <- Norm[j]/decision[, j]
    }
  }
  W <- diag(weights)
  NW <- N %*% W
  posI <- apply(NW, 2, max)
  negI <- apply(NW, 2, min)
  distance = function(x, y) {
    sqrt(sum((x - y)^2))
  }
  posDis <- apply(NW, 1, distance, posI)
  negDis <- apply(NW, 1, distance, negI)
  R <- negDis/(negDis + posDis)
  return(data.frame(Alternatives = 1:nrow(decision), R = R, 
                    Ranking = rank(-R, ties.method = "first")))
}


# VIKOR ----
VIKOR_2 <- function (decision, weights, cb, v) 
{
  if (!is.matrix(decision)) 
    stop("'decision' must be a matrix with the values of the alternatives")
  if (missing(weights)) 
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if (!all.equal(sum(weights),1)) 
    stop("The sum of 'weights' is not equal to 1")
  if (!is.character(cb)) 
    stop("'cb' must be a character vector with the type of the criteria")
  if (!all(cb == "max" | cb == "min")) 
    stop("'cb' should contain only 'max' or 'min'")
  if (length(weights) != ncol(decision)) 
    stop("length of 'weights' does not match the number of the criteria")
  if (length(cb) != ncol(decision)) 
    stop("length of 'cb' does not match the number of the criteria")
  if (missing(v)) 
    stop("a value for 'v' in [0,1] should be provided")
  posI <- as.integer(cb == "max") * apply(decision, 2, max) + 
    as.integer(cb == "min") * apply(decision, 2, min)
  negI <- as.integer(cb == "min") * apply(decision, 2, max) + 
    as.integer(cb == "max") * apply(decision, 2, min)
  norm = function(x, w, p, n) {
    w * ((p - x)/(p - n))
  }
  SAux <- apply(decision, 1, norm, weights, posI, negI)
  S <- apply(SAux, 2, sum)
  R <- apply(SAux, 2, max)
  if (v == 0) 
    Q <- (R - min(R))/(max(R) - min(R))
  else if (v == 1) 
    Q <- (S - min(S))/(max(S) - min(S))
  else Q <- v * (S - min(S))/(max(S) - min(S)) + (1 - v) * 
    (R - min(R))/(max(R) - min(R))
  if (any(Q == "NaN") | any(Q == "Inf")) {
    RankingQ <- rep("-", nrow(decision))
  }
  else {
    RankingQ <- rank(Q, ties.method = "first")
  }
  return(data.frame(Alternatives = 1:nrow(decision), S = S, 
                    R = R, Q = Q, Ranking = RankingQ))
}

# MMOORA ----
MMOORA_2 <- function (decision, weights, cb) 
{
  if (!is.matrix(decision)) 
    stop("'decision' must be a matrix with the values of the alternatives")
  if (missing(weights)) 
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if (!all.equal(sum(weights),1)) 
    stop("The sum of 'weights' is not equal to 1")
  if (!is.character(cb)) 
    stop("'cb' must be a character vector with the type of the criteria")
  if (!all(cb == "max" | cb == "min")) 
    stop("'cb' should contain only 'max' or 'min'")
  if (length(weights) != ncol(decision)) 
    stop("length of 'weights' does not match the number of the criteria")
  if (length(cb) != ncol(decision)) 
    stop("length of 'cb' does not match the number of the criteria")
  d = sqrt(colSums(decision^2))
  NW <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    NW[, j] <- (decision[, j]/d[j]) * weights[j]
  }
  NR <- NW
  for (j in 1:ncol(decision)) {
    if (cb[j] == "min") {
      NR[, j] <- NW[, j] * (-1)
    }
  }
  RS <- apply(NR, 1, sum)
  Ref <- as.integer(cb == "max") * apply(NW, 2, max) + as.integer(cb == 
                                                                    "min") * apply(NW, 2, min)
  RefP <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    RefP[, j] <- abs(Ref[j] - NW[, j])
  }
  RP <- apply(RefP, 1, max)
  max <- NW
  min <- NW
  for (j in 1:ncol(NW)) {
    if (cb[j] == "max") {
      min[, j] <- 1
    }
    else {
      max[, j] <- 1
    }
  }
  A <- apply(max, 1, prod)
  B <- apply(min, 1, prod)
  M <- A/B
  Rrs <- rank(-RS, ties.method = "first")
  Rrp <- rank(RP, ties.method = "first")
  Rm <- rank(-M, ties.method = "first")
  MMRanking = MCDM:::TheoryOfDominance(Rrs, Rrp, Rm, decision)
  return(data.frame(Alternatives = 1:nrow(decision), RatioSystem = RS, 
                    Ranking = Rrs, ReferencePoint = RP, Ranking = Rrp, MultiplicativeForm = M, 
                    Ranking = Rm, MultiMooraRanking = MMRanking))
}

# WASPAS ----
WASPAS_2 <- function (decision, weights, cb, lambda) 
{
  if (!is.matrix(decision)) 
    stop("'decision' must be a matrix with the values of the alternatives")
  if (missing(weights)) 
    stop("a vector containing n weigths, adding up to 1, should be provided")
  if (!all.equal(sum(weights),1))
    stop("The sum of 'weights' is not equal to 1")
  if (!is.character(cb)) 
    stop("'cb' must be a character vector with the type of the criteria")
  if (!all(cb == "max" | cb == "min")) 
    stop("'cb' should contain only 'max' or 'min'")
  if (length(weights) != ncol(decision)) 
    stop("length of 'weights' does not match the number of the criteria")
  if (length(cb) != ncol(decision)) 
    stop("length of 'cb' does not match the number of the criteria")
  if (missing(lambda)) 
    stop("a value for 'lambda' in [0,1] should be provided")
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  Norm <- as.integer(cb == "max") * apply(decision, 2, max) + 
    as.integer(cb == "min") * apply(decision, 2, min)
  N <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    if (cb[j] == "max") {
      N[, j] <- decision[, j]/Norm[j]
    }
    else {
      N[, j] <- Norm[j]/decision[, j]
    }
  }
  W <- diag(weights)
  NW <- N %*% W
  WSM <- apply(NW, 1, sum)
  WPM <- matrix(nrow = nrow(decision), ncol = ncol(decision))
  for (j in 1:ncol(decision)) {
    WPM[, j] <- N[, j]^weights[j]
  }
  WPM <- apply(WPM, 1, prod)
  Q <- (WSM * lambda) + ((1 - lambda) * WPM)
  return(data.frame(Alternatives = 1:nrow(decision), WSM = WSM, 
                    WPM = WPM, Q = Q, Ranking = rank(-Q, ties.method = "first")))
}


