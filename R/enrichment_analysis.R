#' Prognostic Enrichment with Real Data
#'
#' Evaluating biomarkers for prognostic enrichment of clinical trials using real data
#'
#' @param formula Object of class "formula", in the form "outcome ~ predictors", where the outcome is a binary indicator with a value of 1 in cases and a value of 0 in controls.
#' @param data Data frame containing the outcome and predictors specified in the ``formula'' argument. Observations with a missing value of the outcome or of any predictor are dropped. 
#' @param family Character object or call to the family() function specifying the link function that is passed to 'glm' to estimate a risk score when more than one predictor is specified. Defaults to binomial(link = "logit"), which yields logistic regression.
#' @param reduction.under.treatment Number between 0 and 1 indicating the percent reduction in event rate under treatment that the trial should be able to detect with the specified power
#' @param cost.screening Number indicating the cost of screening a patient to determine trial eligibility, This argument is optional; if both the ``cost.screening'' and ``cost.keeping'' arguments are specified, then the total cost of the trial based on each screening threshold is estimated and returned.
#' @param cost.keeping Number indicating the cost of retaining a patient in the trial after enrolling. This argument is optional; if both the ``cost.screening'' and ``cost.keeping'' arguments are specified, then the total cost of the trial based on each screening threshold is estimated and returned.
#' @param do.bootstrap Logical indicating whether bootstrap 95\% confidence intervals should be computed. Defaults to TRUE.
#' @param n.bootstrap Number of bootstrap samples to draw, if ``do.bootstrap'' is set to TRUE. Defaults to 1000.
#' @param power Number between 0 and 1 giving the power the trial should have to reject the null hypothesis that there is no treatment effect. Defaults to 0.9.
#' @param smooth.roc Logical indicating the ``smooth'' argument passed to the roc() function from the `pROC' package when a single biomarker is given. Defaults to FALSE.
#' @param alpha Number between 0 and 1 giving the type I error rate for testing the null hypothesis that there is no treatment effect.  Defaults to 0.025.
#' @param alternative Character specifying whether the alternative hypothesis is one-sided with a higher event rate in the treatment group (``one.sided'') or two-sided (``two.sided''). Defaults to ``one.sided''.
#' @param selected.biomarker.quantiles Numeric vector specifying the quantiles of the biomarker measured in controls that will be used to screen trial participants. Defaults to 0, 0.05, ..., 0.95. All entries must be between at least 0 and less than 1.
#' @return A list with components
#' \itemize{
#'   \item estimates: A data frame with the following summary measures for each biomarker threshold that is used to screen trial participants: `selected.biomarker.quantiles': quantiles of observed biomarker values used for screening. `biomarker.screening.thresholds': the values of the biomarker corresponding to the quantiles, `event.rate': post-screening event rate, `NNS': The estimated number of patients needed to screen to identify one patient eligible for the trial, `SS': The sample size in a clinical trial enrolling only patients whose biomarker-based disease risk is above the level used for screening, `N.screen': The total number of individuals whose biomarker values are screened to determine whether they should be enrolled in the trial, `N.screen.increase.percentage': Percentage in N.screen relative to a trail that does not based on the biomarker. `total.cost': The estimated total cost of running the trial if the biomarker were used for prognostic enrichment (if cost.screening and cost.keeping are specified), `cost.reduction.percentage': The reduction in total cost relative to a trial that does not screen based on the biomarker. 
#'   \item estimates.min.total.cost: The row of the estimates data frame corresponding the screening strategy that results in the lowest total cost.
#'   \item bootstrap.CIs: 95\% bootstrap-based CIs for reported summary measures (if do.bootstrap=TRUE).
#'   \item simulation: A logical indicating whether data were simulated.
#'   \item bootstrap.CIs: 95\% bootstrap-based CIs for reported summary measures (if do.bootstrap=TRUE).
#'   \item biomarker: Biomarker from the given dataset, either the single biomarker specified or the predicted values from logistic regression if multiple biomarkers are specified).
#'   \item response: Response variable specified in the dataset.
#' }
#'
#' @seealso \code{\link{enrichment_simulation}}, \code{\link{plot_enrichment_summaries}}
#' @examples
#' data(dcaData)
#'
#' ## using a single biomarker in the dataset
#' analysis.single.marker <- enrichment_analysis(Cancer ~ Marker1,
#' data=dcaData,
#' reduction.under.treatment=0.3,
#' cost.screening=100,
#' cost.keeping=1000)
#' head(analysis.single.marker$estimates)
#' head(analysis.single.marker$bootstrap.CIs)
#' 
#' ## combining two biomarkers in the dataset
#' analysis.two.markers <- enrichment_analysis(Cancer ~ Marker1 + Marker2,
#' data=dcaData,
#' reduction.under.treatment=0.3,
#' cost.screening=100,
#' cost.keeping=1000)
#' head(analysis.two.markers$estimates)
#' head(analysis.two.markers$bootstrap.CIs)
#' @export
enrichment_analysis <- function(formula,
                                   data,
                                   family=binomial(link=logit),
                                   reduction.under.treatment,
                                   cost.screening=NULL,
                                   cost.keeping=NULL,
                                   do.bootstrap=TRUE,
                                   n.bootstrap=1000,
                                   smooth.roc=FALSE, 
                                   power=0.9,
                                   alpha=0.025,
                                   alternative=c("one.sided", "two.sided"),
                                   selected.biomarker.quantiles=seq(from=0, to=0.95, by=0.05)) {
    alternative <- match.arg(alternative)
    stopifnot(class(formula) == "formula")
    if (!is.data.frame(data)) {
        stop("dataset must be a data frame")
    }
    formula.vars <- all.vars(formula)
    count.vars <- length(formula.vars)
    response.name <- formula.vars[1]
    response <- data[, response.name]
    features.names <- formula.vars[2:count.vars]
    features <- as.matrix(data[, features.names])
    ind.missing <- apply(data[, formula.vars], 1, function(x) sum(is.na(x)) > 0)
    count.any.missing <- sum(ind.missing)
    if (count.any.missing > 0) {
        data <- data[!ind.missing, ]
        warning(paste(count.any.missing, "observation(s) with missing data were removed"))
    }
    if (!all(formula.vars %in% names(data))) {
        stop(paste("Variables named", formula.vars[which(formula.vars %in% names(data) == FALSE)], "were not found in your data"))
    }
    if (!all(response == 0 | response == 1)) {
        stop("The response variable should be binary")
    }
    baseline.event.rate <- mean(response)
    if (!(power > 0 & power < 1)) {
        stop("power should be between 0 and 1")
    }
    if (!(alpha > 0 & alpha < 1)) {
        stop("alpha should be between 0 and 1")
    }
    if (!(reduction.under.treatment > 0 & reduction.under.treatment < 1)) {
        stop("reduction.under.treatment should be between 0 and 1")
    }
    if (!all(selected.biomarker.quantiles >= 0 & selected.biomarker.quantiles < 1)) {
        stop("quantiles of the biomarker measured in controls must be at least 0 and less than 1")
    }
    ##########################################################
    ##  If a dataset is provided, calculate ROC directly (1 biomarker) or estimate with logistic model (> 1 biomarker) ##
    ##########################################################
    if (!is.null(data)) {
        if (count.vars == 2) {
            biomarker.name <- formula.vars[2]
            biomarker <- data[, biomarker.name]
            my.roc <- do.call(what=roc, args=list("formula"=formula, "data"=data, "smooth"=smooth.roc))          
        } else if (count.vars > 2) {
            glm.multiple.markers <- do.call(what=glm, args=list("formula"=formula, "data"=data, "family"=family))
            biomarker <- predict(glm.multiple.markers, type="response")
        }
    }
    biomarker.screening.thresholds <- quantile(biomarker, prob=selected.biomarker.quantiles)
    N <- nrow(data)
    NNS <- sapply(biomarker.screening.thresholds, function(x) N / sum(biomarker > x))
    event.rate <- sapply(biomarker.screening.thresholds, function(x) sum(response[biomarker > x]) / sum(biomarker > x))
    # check for NaN/NA values in event rate (indicates no events for observations above a certain level of biomarker)
    idx.missing.event.rate <- is.na(event.rate)
    if (any(idx.missing.event.rate) == TRUE) {
        warning(paste("There were no events for observations with biomarker values above the following quantile(s): ", do.call(paste, c(as.list(selected.biomarker.quantiles[idx.missing.event.rate]), sep=", ")), ". As a result, we removed ", length(sum(idx.missing.event.rate)), " quantile(s) from the analysis.", sep=""))
        selected.biomarker.quantiles <- selected.biomarker.quantiles[-which(idx.missing.event.rate)]
        biomarker.screening.thresholds <- quantile(biomarker, prob=selected.biomarker.quantiles)
        NNS <- sapply(biomarker.screening.thresholds, function(x) N / sum(biomarker > x))
        event.rate <- sapply(biomarker.screening.thresholds, function(x) sum(response[biomarker > x]) / sum(biomarker > x))
    }
    SS <- sample_size(event.rate=event.rate, reduction.under.treatment=reduction.under.treatment, alpha=alpha, power=power)
    N.screen <- SS * NNS 
    N.screen.increase.percentage <- ((N.screen - N.screen[1]) / N.screen[1]) * 100
    if (!is.null(cost.screening) & !is.null(cost.keeping)) {
        total.cost <- SS * (cost.keeping + cost.screening * NNS)  
        total.cost[1] <- SS[1] * cost.keeping 
        cost.reduction.percentage <- ((total.cost[1] - total.cost) / total.cost[1]) * 100
        total.cost.no.screening <- cost.keeping * SS[1] 
        ind.min.total.cost <- which.min(total.cost)
    } 
    # also allow for bootstrap to estimate standard errors
    if (do.bootstrap == TRUE) {
        n.quantiles <- length(biomarker.screening.thresholds)
        NNS.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        event.rate.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        SS.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        N.screen.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        N.screen.increase.percentage.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        total.cost.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        cost.reduction.percentage.boot <- matrix(NA, nrow=n.quantiles, ncol=n.bootstrap)
        total.cost.no.screening.boot <- rep(NA, n.bootstrap)
        ind.min.total.cost.boot <- rep(NA, n.bootstrap)
        zero.event.rate.count <- 0
        for (b in 1:n.bootstrap) {
            zero.event.rates <- 100 ## just initializing
            while (zero.event.rates > 0) {
                idx.bootstrap <- sample(1:N, size=N, replace=TRUE)
                data.boot <- data[idx.bootstrap, ]
                response.boot <- data.boot[, response.name]
                if (count.vars == 2) {
                    biomarker.boot <- data.boot[, biomarker.name]
                } else if (count.vars > 2) {
                    biomarker.boot <- predict(glm.multiple.markers, newdata=data.boot, type="link")
                }
                biomarker.screening.thresholds.boot <- quantile(biomarker.boot, prob=selected.biomarker.quantiles)
                event.rate.boot[, b] <- sapply(biomarker.screening.thresholds.boot, function(x) sum(response.boot[biomarker.boot > x]) / sum(biomarker.boot > x))
                zero.event.rates <- sum(event.rate.boot[, b] == 0 | is.na(event.rate.boot[, b]))
                if (zero.event.rates > 0) {
                    zero.event.rate.count <- zero.event.rate.count + 1
                }
            }      
            NNS.boot[, b] <- sapply(biomarker.screening.thresholds.boot, function(x) N / sum(biomarker.boot > x))
            SS.boot[, b] <- sample_size(event.rate=event.rate.boot[, b], reduction.under.treatment=reduction.under.treatment, alpha=alpha, power=power)
            N.screen.boot[, b] <- SS.boot[, b] * NNS.boot[, b] # total number of patients needed to be screened
            N.screen.increase.percentage.boot[, b] <- ((N.screen.boot[, b] - N.screen.boot[1, b]) / N.screen.boot[1, b]) * 100
            if (!is.null(cost.screening) & !is.null(cost.keeping)) {
                total.cost.boot[, b] <- SS.boot[, b] * (cost.keeping + cost.screening * NNS.boot[, b]) # total cost
                total.cost.boot[1, b] <- SS.boot[1, b] * cost.keeping # no screening
                cost.reduction.percentage.boot[ , b] <- ((total.cost.boot[1, b] - total.cost.boot[, b]) / total.cost.boot[1, b]) * 100
                total.cost.no.screening.boot[b] <- cost.keeping * SS.boot[1, b]
                ind.min.total.cost.boot[b] <- which.min(total.cost.boot[, b])
            }
        }
        ## get bootstrap CIs
        boot.ci.SS <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.SS) <- c("SS.LB", "SS.UB")
        boot.ci.NNS <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.NNS) <- c("NNS.LB", "NNS.UB")
        boot.ci.event.rate <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.event.rate) <- c("event.rate.LB", "event.rate.UB")
        boot.ci.N.screen <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.N.screen) <- c("N.screen.LB", "N.screen.UB")
        boot.ci.N.screen.increase.percentage <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.N.screen.increase.percentage) <- c("N.screen.increase.percentage.LB", "N.screen.increase.percentage.UB")
        boot.ci.total.cost <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.total.cost) <- c("total.cost.LB", "total.cost.UB")
        boot.ci.cost.reduction.percentage <- matrix(NA, nrow=n.quantiles, ncol=2)
        colnames(boot.ci.cost.reduction.percentage) <- c("cost.reduction.percentage.LB", "cost.reduction.percentage.UB")
        for (r in 1:n.quantiles) {
            boot.ci.SS[r, ] <- quantile(SS.boot[r, ], probs=c(0.025, 0.975))
            boot.ci.NNS[r, ] <- quantile(NNS.boot[r, ], probs=c(0.025, 0.975))
            boot.ci.event.rate[r, ] <- quantile(event.rate.boot[r, ], probs=c(0.025, 0.975))
            boot.ci.N.screen[r, ] <- quantile(N.screen.boot[r, ], probs=c(0.025, 0.975))
            boot.ci.N.screen.increase.percentage[r, ] <- quantile(N.screen.increase.percentage.boot[r, ], probs=c(0.025, 0.975))
            if (!is.null(cost.screening) & !is.null(cost.keeping)) {
                boot.ci.total.cost[r, ] <- quantile(total.cost.boot[r, ], probs=c(0.025, 0.975))
                boot.ci.cost.reduction.percentage[r, ] <- quantile(cost.reduction.percentage.boot[r, ], probs=c(0.025, 0.975))
            }
        }
        if (zero.event.rate.count > 0) {
            warning(paste(zero.event.rate.count, "bootstrap replications had zero events above a specified biomarker threshold and were re-sampled."))
        }
        if (!is.null(cost.screening) & !is.null(cost.keeping)) {
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage, total.cost, cost.reduction.percentage), row.names=NULL)
            estimates$selected.biomarker.quantiles <- estimates$selected.biomarker.quantiles * 100
            rownames(estimates) <- NULL
            bootstrap.CIs <- as.data.frame(cbind(boot.ci.event.rate, boot.ci.NNS, boot.ci.SS, boot.ci.N.screen, boot.ci.N.screen, boot.ci.N.screen.increase.percentage, boot.ci.total.cost, boot.ci.cost.reduction.percentage))
            return(list("biomarker"=biomarker, "response"=response, "simulation"=FALSE, "estimates"=estimates, "bootstrap.CIs"=bootstrap.CIs, "estimates.min.total.cost"=estimates[ind.min.total.cost, ]))
        } else {
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage), row.names=NULL)
            estimates$selected.biomarker.quantiles <- estimates$selected.biomarker.quantiles * 100
            bootstrap.CIs <- as.data.frame(cbind(boot.ci.event.rate, boot.ci.NNS, boot.ci.SS, boot.ci.N.screen, boot.ci.N.screen.increase.percentage))
            return(list("biomarker"=biomarker, "response"=response, "simulation"=FALSE, "estimates"=estimates, "bootstrap.CIs"=bootstrap.CIs))
        }
   } else {
        if (!is.null(cost.screening) & !is.null(cost.keeping)) {
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage, total.cost, cost.reduction.percentage), row.names=NULL)
            estimates$selected.biomarker.quantiles <- estimates$selected.biomarker.quantiles * 100
            return(list("biomarker"=biomarker, "response"=response, "simulation"=FALSE, "estimates"=estimates, "estimates.min.total.cost"=estimates[ind.min.total.cost, ]))
        } else {
            estimates <- as.data.frame(cbind(selected.biomarker.quantiles, biomarker.screening.thresholds, event.rate, NNS, SS, N.screen, N.screen.increase.percentage), row.names=NULL)
            estimates$selected.biomarker.quantiles <- estimates$selected.biomarker.quantiles * 100
            return(list("biomarker"=biomarker, "response"=response, "simulation"=FALSE, "estimates"=estimates))
        }
   }
}

