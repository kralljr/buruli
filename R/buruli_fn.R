#' Seasonal power
#'
#' \code{get_power} Obtains power for seasonal associations using cosignor (Barnett and Dobson 2010)
#'
#' @param N Number of simulations
#' @param a amplitude
#' @param w frequency
#' @param P phase
#' @param t time
#' @param sd standard deviation of random error
#' @param sig significance level (defaults to 0.05)
#' @export
#' @examples
#' N <- 100
#' a <- 1
#' w <- 1/16
#' P <- 0
#' t <- seq(1, 365 * 2, by = 1)
#' get_power(100, a, w, P, t)
get_power <- function(N, a, w, P, t, sd1 = 1, sig = 0.05) {
    
    #get y
    y <- a * cos(w * t - P) + rnorm(t, sd = sd1)
    
    #get wt
    t1 <- t * (t <= 365) + (t - 365) * (t > 365)
    ft <- (t1 - 1)/365
    wt <- 2 * pi * ft
    
    sig <- vector(, length = N)
    for(i in 1 : N) {
        sig[i] <- get_sig(y, wt, sig)
    }
    
    return(mean(sig, na.rm = T))
}

#' Seasonal power
#'
#' \code{get_sig} Get significance info for one simulation
#'
#' @param y Response
#' @param wt time in cycle
#' @param sig significance level (defaults to 0.05)
get_sig <- function(y, wt, sig = 0.05) {
    glm1 <- glm(y ~ cos(wt) + sin(wt))
    sglm1 <- summary(glm1)$coef[-1, 4]
    
    out <- 0
    if(sglm1[1] < sig / 2 | sglm1[2] < sig / 2) {
        out <- 1
    }
    
    return(out)
}