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
#' @param max1 Maximum number of days delayed
#' @export
#' @examples
#' N <- 100
#' a <- 1
#' w <- 1/16
#' P <- 0
#' t <- seq(1, 365 * 2, by = 1)
#' get_power(100, a, w, P, t)
get_power <- function(N, a, w, P, t, sd = 1, sig = 0.05, max1 = 0) {
    
    #get wt
    t1 <- t * (t <= 365) + (t - 365) * (t > 365)
    ft <- (t1 - 1)/365
    wt <- 2 * pi * ft
    
    sig.out <- vector(, length = N)
    for(i in 1 : N) {
        #get y
        mu <- (  10 + a * cos(wt) + 0 * sin(wt) + rnorm(length(wt), sd = sd)) #/ 10000
        
        #get true data
        #y <- rbinom(length(mu), 10000, mu)
        y <- rpois(length(mu), mu) 
        
        if(max1 > 0) {
            #get individual data info
            dates <- as.Date(t, origin = "1970-01-01")
            
#             y2 <- rep(0, length(y))
#             for(j in 1 : length(y)) {
#                 del <- sample(seq(0, max1), length(y[j]), replace = T)
#                 y2[] <- y2[] + del
#                 for(k in 1 : y[j]) {
#                     
#                 }
#             }
#             
            ind <- as.Date("1970-01-01", origin = "1970-01-01")
            for(j in 1 : length(y)) {
                ind <- c(ind, rep(dates[j], y[j]))
            }
            ind <- ind[-1]
            
            
            #add delay and get count/day
            del <- sample(seq(0, max1), length(ind), replace = T)
            ind <- ind + del
            ind <- dplyr::filter(data.frame(ind), ind <= max(dates, na.rm = T))
            td <- table(as.matrix(ind))
            y <- td[order(as.Date(names(td), origin = "1970-01-01"))]
            
            #add in missing dates
            dates <- data.frame(dates)
            y <- data.frame(as.Date(names(y), origin = "1970-01-01"), y)
            colnames(y) <- c("dates", "y")
            y <- merge(y, dates, all.y = T)
            y[is.na(y)] <- 0
            y <- y$y
        }

        sig.out[i] <- get_sig(y, wt, sig)
    }
    
    out <- mean(sig.out, na.rm = T)
    out
    return(out)
}

#' Seasonal power
#'
#' \code{get_sig} Get significance info for one simulation
#'
#' @param y Response
#' @param wt time in cycle
#' @param sig significance level (defaults to 0.05)
get_sig <- function(y, wt, sig = 0.05) {
    glm1 <- glm(y ~ cos(wt) + sin(wt), family= "poisson")
    sglm1 <- summary(glm1)$coef[-1, 4]
    
    out <- 0
    if(sglm1[1] < sig / 2 | sglm1[2] < sig / 2) {
        out <- 1
    }
    
    return(out)
}