#' Prognostic Enrichment with Simulated Data
#'
#' Evaluating biomarkers for prognostic enrichment of clinical trials using simulated data
#'
#' @param baseline.event.rate  A number between 0 and 1 indicating the prevalence of the event in the study population.
#' @param reduction.under.treatment A number between 0 and 1 indicating the percent reduction in event rate under treatment that the trial should be able to detect with the specified power.
#' @param estimated.auc A numeric vector, with each entry between 0.5 and 1, that specifies the AUC for each biomarker to use in simulations.
#' @param cost.screening A positive number indicating the cost of screening a patient to determine trial eligibility, This argument is optional; if both cost.screening and cost.keeping are specified, then then the total cost of the trial based on each screening threshold is estimated and returned.
#' @param cost.keeping A positive number indicating the cost of retaining a patient in the trial after enrolling. This argument is optional; if both cost.screening and cost.keeping are specified, then then the total cost of the trial based on each screening threshold is estimated and returned.
#' @param roc.type A character vector with the same length as the estimated.auc argument. Each entry must be one of "symmetric", "right.shifted", or "left.shifted", which describes the general shape of the ROC curve to use for simulated data. Defaults to "symmetric" for each biomarker.
#' @param simulation.sample.size A positive number giving the sample size to use for simulated data. Defaults to 500,000 (to help see trends).
#' @param power Number between 0 and 1 giving the power the trial should have to reject the null hypothesis that there is no treatment effect. Defaults to 0.9.
#' @param alpha Number between 0 and 1 giving the type I error rate for testing the null hypothesis that there is no treatment effect.  Defaults to 0.025.
#' @param alternative Character specifying whether the alternative hypothesis is one-sided (``one.sided'') with a higher outcome probability in the treatment group or two-sided (``two.sided''). Defaults to ``one.sided''.
#' @param selected.biomarker.quantiles Numeric vector specifying the quantiles of the biomarker measured in controls that will be used to screen trial participants. Defaults to 0, 5, ..., 95. All entries must be between at least 0 and less than 001.
#' @return A list with components
#' \itemize{
#'   \item estimates: A data frame with the following summary measures for each biomarker threshold that is used to screen trial participants: `selected.biomarker.quantiles': quantiles of observed biomarker values used for screening. `biomarker.screening.thresholds': the values of the biomarker corresponding to the quantiles, `event.rate': post-screening event rate, `NNS': The estimated number of patients needed to screen to identify one patient eligible for the trial, `SS': The sample size in a clinical trial enrolling only patients whose biomarker-based disease risk is above the level used for screening, `N.screen': The total number of individuals whose biomarker values are screened to determine whether they should be enrolled in the trial, `N.screen.increase.percentage': Percentage in N.screen relative to a trail that does not based on the biomarker. `total.cost': The estimated total cost of running the trial if the biomarker were used for prognostic enrichment (if cost.screening and cost.keeping are specified), `cost.reduction.percentage': The reduction in total cost relative to a trial that does not screen based on the biomarker. `Biomarker': label for the biomarker. 
#`   \item roc.data: Data frame with three columns -- 'FPR', 'TPR', and `Biomarker' -- that are used to make the ROC plots for each biomarker.
#'   \item simulation: Logical indicating whether data were simulated (always TRUE for the \code{\link{plot_enrichment_summaries}} function).
#' }
#'
#' @seealso \code{\link{enrichment_analysis}}, \code{\link{plot_enrichment_summaries}}
#' @examples
#' ## three biomarkers with symmetric ROC curves
#' simulation.three.markers <- enrichment_simulation(baseline.event.rate=0.2,
#' reduction.under.treatment=0.3,
#' estimated.auc=c(0.72, 0.82, 0.85),
#' roc.type=c("symmetric", "symmetric", "symmetric"),
#' cost.screening=1,
#' cost.keeping=10,
#' simulation.sample.size=1e+5)
#' head(simulation.three.markers$estimates)
#'
#' @import pROC
#' @import VGAM
#' @importFrom graphics abline lines par
#' @importFrom stats binomial glm qnorm quantile rbinom rnorm
#' @export
enrichment_simulation <- function(baseline.event.rate,
                                   reduction.under.treatment,
                                   estimated.auc,
                                   roc.type=NULL,
                                   cost.screening=NULL,
                                   cost.keeping=NULL,
                                   simulation.sample.size=5e+5,
                                   alternative=c("one.sided", "two.sided"),
                                   power=0.9,
                                   alpha=0.025,
                                   selected.biomarker.quantiles=seq(from=0, to=95, by=5)) {

    #############
    ## Check arguments ##
    #############
    stopifnot(is.numeric(baseline.event.rate))
    stopifnot(baseline.event.rate > 0 & baseline.event.rate < 1)
    updated.auc <- estimated.auc[!is.na(estimated.auc)]
    n.auc <- length(updated.auc)
    for (i in 1:n.auc) {
        stopifnot(updated.auc[i] >= 0.5 & updated.auc[i] <= 1)
    }
    if (is.null(roc.type)) {
        roc.type <- rep("symmetric", n.auc)
    }
    n.roc.type <- length(roc.type)
    if (n.auc != n.roc.type) {
        stop("estimated.auc and roc.type need to have the same length")
    }
    if (is.null(cost.screening) & !is.null(cost.keeping)) {
        stop("if cost.screening is specified, then cost.keeping must also be specified")
    }
    if (!is.null(cost.screening) & is.null(cost.keeping)) {
        stop("if cost.keeping is specified, then cost.screening must also be specified")
    }
    if (!is.null(cost.screening) & !is.null(cost.keeping) & (cost.screening < 0 | cost.keeping < 0)) {
        stop("cost.keeping and cost.keeping should both be positive numbers")
    }
    for (i in 1:n.roc.type) {
        if (!(roc.type[i] %in% c("symmetric", "right.shifted", "left.shifted"))) {
            stop("each entry of the roc.type argument needs to be either symmetric, right.shifted, or left.shifted")
        }
    }
    updated.roc.type <- roc.type[!is.na(estimated.auc)]
    alternative <- match.arg(alternative)
    stopifnot(is.numeric(power))
    if (!(power > 0 & power < 1)) {
        stop("power should be between 0 and 1")
    }
    if (!(alpha > 0 & alpha < 1)) {
        stop("alpha should be between 0 and 1")
    }
    if (!(reduction.under.treatment > 0 & reduction.under.treatment < 1)) {
        stop("reduction.under.treatment should be between 0 and 1")
    }
    if (!all(selected.biomarker.quantiles >= 0 & selected.biomarker.quantiles < 100)) {
        stop("quantiles of the biomarker measured in controls must be at least 0 and less than 100")
    }
    if (simulation.sample.size < 0) {
        stop("simulation.sample.size should be a positive number")
    }
    N <- simulation.sample.size
    estimates.check <- NULL
    roc.data.check <- NULL
    for (i in 1:n.auc) {
        simulation.data <- user_auc_and_roc_type_to_data(N=N, baseline.event.rate=baseline.event.rate,
                                                                                auc=updated.auc[i], roc.type=updated.roc.type[i],
                                                                                selected.biomarker.quantiles=selected.biomarker.quantiles / 100)
        roc.data <- user_auc_to_plots(auc=updated.auc[i], baseline.event.rate=baseline.event.rate, roc.type=updated.roc.type[i], prototypical=FALSE)
        biomarker <- simulation.data$biomarker
        response <- simulation.data$response
        biomarker.screening.thresholds <- quantile(biomarker, prob=selected.biomarker.quantiles / 100)
        if (updated.roc.type[i] %in% c("symmetric", "left.shifted")) {
            NNS <- sapply(biomarker.screening.thresholds, function(x) N / sum(biomarker > x))
            event.rate <- sapply(biomarker.screening.thresholds, function(x) sum(response[biomarker > x]) / sum(biomarker > x))
        } else {
            NNS <- sapply(biomarker.screening.thresholds, function(x) N / sum(-biomarker > x))
            event.rate <- sapply(biomarker.screening.thresholds, function(x) sum(response[-biomarker > x]) / sum(-biomarker > x))
        }
        SS <- sample_size(event.rate=event.rate, reduction.under.treatment=reduction.under.treatment, alpha=alpha, power=power, alternative=alternative)
        N.screen <- SS * NNS 
        N.screen.increase.percentage <- ((N.screen - N.screen[1]) / N.screen[1]) * 100
        cost.missing <- is.null(cost.screening) | is.null(cost.keeping)
        if (cost.missing == TRUE) {
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage))
            table.data <- as.data.frame(cbind(paste(selected.biomarker.quantiles, "%", sep=""),
                                              round(event.rate, 2),
                                              round(SS, 0),
                                              round(NNS, 1),
                                              round(N.screen, 0),
                                              round(N.screen.increase.percentage, 1)))
            table.data <- as.data.frame(cbind(apply(table.data[, 1:5], 2, function(x) as.character(x)),
                                              paste(as.character(table.data[, 6]), "%", sep="")))
            names(table.data) <- c("Percent of Patients Screened from Trial", "Event Rate Among Biomarker-Positive Patients", "Sample Size", "NNS", "Total Screened", "Percent Change in Total Screened")
            rownames(table.data) <- NULL
        } else {
            total.cost <- SS * (cost.keeping + cost.screening * NNS) # total cost
            total.cost[1] <- SS[1] * cost.keeping
            cost.reduction.percentage <- ((total.cost[1] - total.cost) / total.cost[1]) * 100
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage, total.cost, cost.reduction.percentage))
            table.data <- as.data.frame(cbind(paste(selected.biomarker.quantiles , "%", sep=""),
                                              round(event.rate, 2),
                                              round(SS, 0),
                                              round(NNS, 1),
                                              round(N.screen, 0),
                                              round(N.screen.increase.percentage, 1),
                                              round(total.cost, 0),
                                              round(cost.reduction.percentage, 1)))
            table.data <- cbind(apply(table.data[, 1:5], 2, function(x) as.character(x)),
                                paste(as.character(table.data[, 6]), "%", sep=""),
                                as.character(table.data[, 7]),
                                paste(as.character(table.data[, 8]), "%", sep=""))
            table.data <- as.data.frame(table.data)
            names(table.data) <- c("Percent of Patients Screened from Trial", "Event Rate Among Biomarker-Positive Patients", "Sample Size", "NNS", "Total Screened", "Percent Change in Total Screened", "Total Costs for Screening and Patients in Trial", "Percent Reduction in Total Cost")
            rownames(table.data) <- NULL
        }
        estimates$Biomarker <- paste("Biomarker", as.character(i), sep=" ")
        estimates.check <- rbind(estimates.check, estimates)
        roc.df <- as.data.frame(cbind(roc.data$fpr.vec, roc.data$tpr.vec))
        roc.df$Biomarker <- paste("Biomarker", as.character(i), sep=" ")
        roc.data.check <- rbind(roc.data.check, roc.df)
    }
    names(roc.data.check) <- c("FPR", "TPR", "Biomarker")
    return(list("simulation"=TRUE, "estimates"=estimates.check,"roc.data"=roc.data.check))
}

user_auc_and_roc_type_to_data <- function(N, baseline.event.rate, auc, roc.type, selected.biomarker.quantiles) {
    if (roc.type == "right.shifted") {
        # need to flip case/control labels and eventually call a test "positive" if it is below the threshold, rather than above
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        biomarker <- numeric(N)
        biomarker[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=1)
        biomarker[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=(1-auc)/auc)
        biomarker.screening.thresholds <- quantile(-biomarker, prob=selected.biomarker.quantiles)
    }
    if (roc.type == "symmetric") {
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        sd.cases <- 1
        b <- 1 / sd.cases
        a <- sqrt(1 + b^2) * qnorm(auc)
        mean.cases <- a * sd.cases 
        biomarker <- numeric(N)
        biomarker[response == 1] <- rnorm(n=sum(response==1), mean=mean.cases, sd=sd.cases)
        biomarker[response == 0] <- rnorm(n=sum(response==0), mean=0, sd=1)
        biomarker.screening.thresholds <- quantile(biomarker, prob=selected.biomarker.quantiles)
    } else if (roc.type == "left.shifted") {
        response <- rbinom(n=N, size=1, prob=baseline.event.rate)
        biomarker <- numeric(N)
        biomarker[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=1)
        biomarker[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=(1-auc)/auc)
        biomarker.screening.thresholds <- quantile(biomarker, prob=selected.biomarker.quantiles)
    }
    return(list("selected.biomarker.quantiles"=selected.biomarker.quantiles,
                   "biomarker.screening.thresholds"=biomarker.screening.thresholds,
                   "response"=response, "biomarker"=biomarker))
}

user_auc_to_plots <- function(auc, baseline.event.rate, roc.type=NULL, prototypical=TRUE, n=5e+4, n.thresholds=1000, verbose=TRUE) {
    response <- rbinom(n, size=1, prob=baseline.event.rate)
    ## high TPR earlier (orange)
    x.left.shifted <- rep(NA, n)
    x.left.shifted[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=1)
    x.left.shifted[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=(1-auc)/auc)
    result.left.shifted <- get_roc(x=x.left.shifted, response=response, test.positive="higher", n.thresholds=n.thresholds, verbose=verbose)
    ## high TPR earlier (blue) 
    x.right.shifted <- rep(NA, n)
    x.right.shifted[response == 1] <- rlomax(n=sum(response==1), scale=1, shape3.q=1)
    x.right.shifted[response == 0] <- rlomax(n=sum(response==0), scale=1, shape3.q=(1-auc)/auc)
    result.right.shifted <- get_roc(x=x.right.shifted, response=response, test.positive="lower", n.thresholds=n.thresholds, verbose=verbose)
    ## symmetric (black)
    x.symmetric <- rep(NA, n)
    x.symmetric[response == 0] <- rnorm(n=sum(response==0), mean=0, sd=1)
    x.symmetric[response == 1] <- rnorm(n=sum(response==1), mean=sqrt(2) * qnorm(auc), sd=1)
    result.symmetric <- get_roc(x=x.symmetric, response=response, test.positive="higher", n.thresholds=n.thresholds, verbose=verbose)
    # make plots
    if (prototypical == TRUE) {
        par(pty="s")
        plot(result.left.shifted$fpr.vec, result.left.shifted$tpr.vec, type="l", lty=2, lwd=3, col="orange", ylim=c(0, 1), xlim=c(0, 1),
             xlab="FPR (%)", ylab="TPR (%)",
             main=paste("Prototypical ROC Curves \n (Example for AUC=", auc, ")", sep=""))
        lines(result.right.shifted$fpr.vec, result.right.shifted$tpr.vec, type="l", lty=3, lwd=3, col="cyan")
        lines(result.symmetric$fpr.vec, result.symmetric$tpr.vec, type="l", lty=1, lwd=3, col="black")
        abline(a=0, b=1, lwd=2, lty="dashed", col="gray")
    } else {
        if (roc.type == "symmetric") {
            fpr.vec <- result.symmetric$fpr.vec
            tpr.vec <- result.symmetric$tpr.vec
        } else if (roc.type == "right.shifted") {
            fpr.vec <- result.right.shifted$fpr.vec
            tpr.vec <- result.right.shifted$tpr.vec
        } else if (roc.type == "left.shifted") {
            fpr.vec <- result.left.shifted$fpr.vec
            tpr.vec <- result.left.shifted$tpr.vec
        }
        return(list("fpr.vec"=fpr.vec, "tpr.vec"=tpr.vec))
    }
}


get_roc <- function(x, response, test.positive=c("higher","lower"), n.thresholds=1000, verbose=TRUE) {
    test.positive <- match.arg(test.positive)
    n.healthy <- sum(response == 0)
    n.disease <- sum(response == 1)
    if (test.positive == "lower") {
        x <- -x
    }
    thresholds <- as.numeric(sort(quantile(x=x[response==0], prob= seq(from=0, to=0.99, length.out=n.thresholds), decreasing=FALSE)))
    tpr <- rep(NA, n.thresholds)
    fpr <- rep(NA, n.thresholds)
    x.cases <- x[response==1]
    x.controls <- x[response==0]
    for (i in 1:n.thresholds) {
        tpr[i] <- sum(x.cases > thresholds[i]) / n.disease
        fpr[i] <- sum(x.controls > thresholds[i]) / n.healthy
    }
    if (verbose == TRUE) {
        auc.estimate <- mean(sample(x.cases, size=5e+5, replace=TRUE) > sample(x.controls, size=5e+5, replace=TRUE))
        return(list("auc.estimate"=auc.estimate, "fpr.vec"=fpr, "tpr.vec"=tpr))
    } else {
        return(list("fpr.vec"=fpr, "tpr.vec"=tpr))
    }
}

