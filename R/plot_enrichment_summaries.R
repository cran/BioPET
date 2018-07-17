#' Plotting Prognostic Enrichment Estimates
#'
#' Plot summaries of prognostic enrichment of clinical trials estimated by the \code{\link{enrichment_analysis}} and \code{\link{enrichment_simulation}} functions.
#' @param x Object returned by either the \code{\link{enrichment_analysis}} or the \code{\link{enrichment_simulation}} function.
#' @param text.size.x.axis Size of text for the x-axis of plots. Defaults to 10.
#' @param text.size.y.axis Size of text for the y-axis of plots. Defaults to 10.
#' @param text.size.plot.title Size of text for the plot titles. Defaults to 10.
#' @param text.size.axis.ticks Size of axis tick marks for plots. Defaults to 10.
#' @param annotate.no.screening.cost Logical indicating whether to annotate the relative to total cost curve at the point where no biomarker screening occurs. Defaults to FALSE
#' @param smooth.roc  Logical indicating whether the ROC curves (plotting with the roc() function in the `pROC' package) should be smoothed. Defaults to TRUE.
#' @return A grid of either 4 or 6 plots, summarizing the results of either the \code{\link{enrichment_analysis}} or the \code{\link{enrichment_simulation}} function.
#' @seealso \code{\link{enrichment_analysis}}, \code{\link{enrichment_simulation}}
#' @examples
#' 
#' data(dcaData)
#' # one marker
#' analysis.single.marker <- enrichment_analysis(Cancer ~ Marker1,
#' data=dcaData,
#' reduction.under.treatment=0.3,
#' cost.screening=100, cost.keeping=1000)
#' plot_enrichment_summaries(analysis.single.marker)
#'
#' # two markers
#' analysis.two.markers <- enrichment_analysis(Cancer ~ Marker1 + Marker2,
#' data=dcaData,
#' reduction.under.treatment=0.3,
#' cost.screening=100,
#' cost.keeping=1000)
#' plot_enrichment_summaries(analysis.two.markers)
#' @import ggplot2
#' @import pROC
#' @import gridExtra
#' @export
plot_enrichment_summaries <- function(x,
                                           text.size.x.axis=10,
                                           text.size.y.axis=10,
                                           text.size.plot.title=10,
                                           text.size.axis.ticks=10,
                                           annotate.no.screening.cost=FALSE,
                                           smooth.roc=TRUE) {
    
    ## Determine whether bootstrap CIs are available and store data used for plotting
    bootstrap.indicator <- "bootstrap.CIs" %in% names(x)
    if (bootstrap.indicator == TRUE) {
        estimates <- cbind(x$estimates, x$bootstrap.CIs)
    } else {
        estimates <- x$estimates
    }
    ## Determine whether we are using real or simulated data
    using.simulated.data <- x$simulation
    if (using.simulated.data == TRUE) {
        roc.data <- x$roc.data
        roc.data$FPR <- roc.data$FPR * 100
        roc.data$TPR <- roc.data$TPR * 100
    } else {
        real.roc <- roc(response=x$response, predictor=x$biomarker, smooth=smooth.roc)
        roc.data <- as.data.frame(cbind((1 - real.roc$specificities) * 100, real.roc$sensitivities * 100))
        names(roc.data) <- c("FPR", "TPR")
    }
    ## Determine whether the total cost information is in the data frame (i.e. whether the user specified cost.screening and cost.retention earlier)
    cost.indicator <- "total.cost" %in% names(estimates) & "cost.reduction.percentage" %in% names(estimates)
    ## Plot ROC curves for specified biomarkers (simulated data only)
    if (using.simulated.data == TRUE) {
        plot.ROC <- ggplot(roc.data, aes_string(x="FPR", y="TPR", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) + geom_line(size=0.9) + geom_abline(intercept = 0, slope = 1, linetype="dashed", colour="gray") + coord_fixed() +
            labs(x="FPR (%)", y="TPR (%)") +
            ggtitle("ROC Curve for Specified Biomarkers")
    } else {
        plot.ROC <- ggplot(roc.data, aes_string(x="FPR", y="TPR")) + geom_line(size=0.9) + geom_abline(intercept = 0, slope = 1, linetype="dashed", colour="gray") + coord_fixed() +
            labs(x="FPR (%)", y="TPR (%)") +
            ggtitle("ROC Curve for Specified Biomarker")
    }
    theme_update(plot.title=element_text(hjust=0.5))
    plot.ROC <- plot.ROC + expand_limits(y=0) +
        theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title))  + scale_x_continuous(expand = c(0, 0)) +
        theme(legend.title=element_blank())
    if (cost.indicator == TRUE) {
        plot.ROC <- plot.ROC + theme(legend.text=element_text(size=10), legend.key.size = unit(0.45, "cm"), legend.position=c(0.6, 0.25))
    } else {
        plot.ROC <- plot.ROC + theme(legend.text=element_text(size=13), legend.key.size = unit(0.6, "cm"), legend.position=c(0.7, 0.2))
    }
    ## Biomarker percentile vs. sample size
    if (using.simulated.data == TRUE) {
        plot.sample.size <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="SS", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) +
            labs(x="Percent of Patients Screened from Trial", y="Sample Size") +
            ggtitle("Clinical Trial Total Sample Size")
    } else {
        plot.sample.size <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="SS")) +
            labs(x="Percent of Patients Screened from Trial", y="Sample Size") +
            ggtitle("Clinical Trial Total Sample Size")
    }
    plot.sample.size <- plot.sample.size + geom_line(size=0.9) + geom_point(size=1.5) + theme(axis.title.x = element_text(size=text.size.x.axis)) +
        theme(axis.title.y = element_text(size=text.size.y.axis)) +
        theme(axis.text= element_text(size=text.size.axis.ticks)) +
        theme(plot.title = element_text(size=text.size.plot.title)) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
    if (cost.indicator == TRUE) {
        plot.sample.size <-  plot.sample.size + labs(x="", y="Total Screened")
    } else {
        plot.sample.size <- plot.sample.size + labs(x="Percent of Patients Screened from Trial", y="Total Screened")
    }
    ## Biomarker percentile vs. event rate after screening
    if (using.simulated.data == TRUE) {
        plot.event.rate <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="event.rate", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) +
                          labs(x="", y="Event Rate") +
                           ggtitle("Event Rate Among \n Biomarker-Positive Patients")
    } else {
        plot.event.rate <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="event.rate")) + 
           labs(x="", y="Post-Screening Event Rate") + 
            ggtitle("Event Rate Among \n Biomarker-Positive Patients")
    }
        plot.event.rate <- plot.event.rate + geom_line(size=0.9) + geom_point(size=1.5) + expand_limits(y=0) +
            theme(axis.title.x = element_text(size=text.size.x.axis)) +
            theme(axis.title.y = element_text(size=text.size.y.axis)) +
            theme(axis.text= element_text(size=text.size.axis.ticks)) +
            theme(plot.title = element_text(size=text.size.plot.title)) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
    ## Biomarker percentile vs. total # needing to be screened
    if (using.simulated.data == TRUE) {
        plot.N.screen <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="N.screen", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) +
            ggtitle("Total Number of Patients \n Screened to Enroll Trial") 
    } else {
        plot.N.screen <- ggplot(estimates, aes_string(x="selected.biomarker.quantiles", y="N.screen")) + geom_line() + geom_point(size=1.5) +
            ggtitle("Total Number of Patients \n Screened to Enroll Trial")
    }
    if (cost.indicator == TRUE) {
        plot.N.screen <-  plot.N.screen + labs(x="", y="Total Screened")
    } else {
        plot.N.screen <- plot.N.screen + labs(x="Percent of Patients Screened from Trial", y="Total Screened")
    }
        plot.N.screen <- plot.N.screen + geom_line(size=0.9) + geom_point(size=1.5) + expand_limits(y=0) +
            theme(axis.title.x = element_text(size=text.size.x.axis)) +
            theme(axis.title.y = element_text(size=text.size.y.axis)) +
            theme(axis.text= element_text(size=text.size.axis.ticks)) +
            theme(plot.title = element_text(size=text.size.plot.title)) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
    
    if (bootstrap.indicator == TRUE) {
        ## store bootstrap limits
        limits.sample.size <- aes_string(ymin="SS.LB", ymax="SS.UB")
        limits.event.rate <- aes_string(ymin="event.rate.LB", ymax="event.rate.UB")
        limits.N.screen <- aes_string(ymin="N.screen.LB", ymax="N.screen.UB")
        ## plot bootstrap CIs
        plot.sample.size <- plot.sample.size + geom_errorbar(limits.sample.size, width=5, linetype="dashed", colour="darkblue") + coord_cartesian(ylim=c(0, estimates$SS.UB[1]*1.05))
        plot.event.rate <- plot.event.rate + geom_errorbar(limits.event.rate, width=5, linetype="dashed", colour="darkblue") + coord_cartesian(ylim=c(0, max(estimates$event.rate.UB)*1.05))
        plot.N.screen <- plot.N.screen + geom_errorbar(limits.N.screen, width=5, linetype="dashed", colour="darkblue") + coord_cartesian(ylim=c(0, max(estimates$N.screen)*1.05))
    }
    ## Biomarker percentile vs. total cost
    if (cost.indicator == TRUE) {
        ind.total.cost <- 1:nrow(estimates)
        if (using.simulated.data == TRUE) {
            plot.total.cost <- ggplot(estimates[ind.total.cost, ], aes_string(x="selected.biomarker.quantiles", y="total.cost", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) + geom_line(size=0.9) + geom_point(size=2) + geom_hline(yintercept=0, linetype=2) +
                labs(x="Percent of Patients Screened from Trial", y="Total Cost") +
                ggtitle("Total Costs for Screening \n and Patients in Trial")
        } else {
            plot.total.cost <- ggplot(estimates[ind.total.cost, ], aes_string(x="selected.biomarker.quantiles", y="total.cost")) + geom_line() + geom_hline(yintercept=0, linetype=2) + geom_point(size=1.5) +
                labs(x="Percent of Patients Screened from Trial", y="Total Cost") +
                ggtitle("Total Costs for Screening \n and Patients in Trial")
        }
        plot.total.cost <- plot.total.cost + expand_limits(y=0)+ 
            theme(axis.title.x = element_text(size=text.size.x.axis)) +
            theme(axis.title.y = element_text(size=text.size.y.axis)) +
            theme(axis.text= element_text(size=text.size.axis.ticks)) +
            theme(plot.title = element_text(size=text.size.plot.title)) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
        if (annotate.no.screening.cost == TRUE) {
            plot.total.cost <- plot.total.cost + annotate("text", x=estimates[1, "selected.biomarker.quantiles"], y=0.99 * estimates[1, "total.cost"],
                                                          angle=90, color="blue", size=2, label="no screening") +
                annotate("point", x=estimates[1, "selected.biomarker.quantiles"], y=estimates[1, "total.cost"],
                         size=2, color="blue")
        }
        ## Biomarker percentile vs. percentage reduction in total cost (relative to no screening scenario)
        if (using.simulated.data == TRUE) {
            plot.cost.reduction.percentage <- ggplot(estimates[ind.total.cost, ], aes_string(x="selected.biomarker.quantiles", y="cost.reduction.percentage", group="Biomarker", shape="Biomarker", col="Biomarker", linetype="Biomarker", fill="Biomarker")) + geom_line(size=0.9) + geom_point(size=1.5) + geom_hline(yintercept=0, linetype=2) +
               labs(x="Percent of Patients Screened from Trial", y="% Reduction in Total Cost") +
                ggtitle("Percent Reduction \n in Total Cost")
        } else {
            plot.cost.reduction.percentage <- ggplot(estimates[ind.total.cost, ], aes_string(x="selected.biomarker.quantiles", y="cost.reduction.percentage")) + geom_line() + geom_hline(yintercept=0, linetype=2) + geom_point(size=1.5) +
               labs(x="Percent of Patients Screened from Trial", y="% Reduction in Total Cost") +
                ggtitle("Percent Reduction \n in Total Cost")
        }
        plot.cost.reduction.percentage <- plot.cost.reduction.percentage + expand_limits(y=0) +
            theme(axis.title.x = element_text(size=text.size.x.axis)) +
            theme(axis.title.y = element_text(size=text.size.y.axis)) +
            theme(axis.text= element_text(size=text.size.axis.ticks)) +
            theme(plot.title = element_text(size=text.size.plot.title)) + scale_x_continuous(expand = c(0, 0)) + theme(legend.position = 'none')
        if (bootstrap.indicator == TRUE) {
            ## bootstrap limits
            limits.total.cost <- aes_string(ymin="total.cost.LB", ymax="total.cost.UB")
            limits.cost.reduction.percentage <- aes_string(ymin="cost.reduction.percentage.LB", ymax="cost.reduction.percentage.UB")
            lower.limit.cost.reduction.percentage <- max(min(estimates$cost.reduction.percentage.LB, 0), -20)
            ## plotting bootstrap CIs
            plot.total.cost <- plot.total.cost + geom_errorbar(limits.total.cost, width=5, linetype="dashed", colour="darkblue") + coord_cartesian(ylim=c(0, estimates$total.cost.UB[1]*1.05))
            plot.cost.reduction.percentage <- plot.cost.reduction.percentage + geom_errorbar(limits.cost.reduction.percentage, width=5, linetype="dashed", colour="darkblue") + coord_cartesian(ylim=c(lower.limit.cost.reduction.percentage, max(estimates$cost.reduction.percentage.UB) * 1.05))
        }
        args.to.plot <- list("plot.ROC"=plot.ROC, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size, "plot.N.screen"=plot.N.screen, "plot.total.cost"=plot.total.cost, "plot.cost.reduction.percentage"=plot.cost.reduction.percentage)
    } else {
        args.to.plot <- list("plot.ROC"=plot.ROC, "plot.event.rate"=plot.event.rate, "plot.sample.size"=plot.sample.size,  "plot.N.screen"=plot.N.screen)
    }
    ## Show the plots the user wants to see in a grid
    do.call(grid.arrange, c(args.to.plot, list(ncol=2)))
}

