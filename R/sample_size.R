sample_size <- function(event.rate, reduction.under.treatment,
                                 alpha=0.025, power=0.9,
                                 alternative=c("one.sided", "two.sided")) {
    alternative <- match.arg(alternative)
    p1 <- event.rate
    p2 <- event.rate * (1 - reduction.under.treatment)
    z1.one.sided <- qnorm(1 - alpha)
    z1.two.sided <- qnorm(1 - alpha/2)
    z2 <- qnorm(power)
    if (alternative == "one.sided") { 
        SS <- 2 * ( (z1.one.sided * sqrt((p1 + p2) * (1 - (p1 + p2)/2) ) +
                        z2 * sqrt(p1 * (1 - p1) + p2 * (1 - p2)) )^2 / (p1 - p2)^2)
    } else if (alternative == "two.sided") {
        SS <- 2 * ( (z1.two.sided * sqrt((p1 + p2) * (1 - (p1 + p2)/2) ) +
                        z2 * sqrt(p1 * (1 - p1) + p2 * (1 - p2)) )^2 / (p1 - p2)^2)
    } 
    return(SS)
}

