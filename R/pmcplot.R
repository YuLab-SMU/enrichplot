##' PubMed Central Trend plot
##'
##'
##' @title pmcplot
##' @param query query terms
##' @param period period of query in the unit of year
##' @param proportion If TRUE, use query_hits/all_hits, otherwise use query_hits
##' @return ggplot object
##' @importFrom purrr map_df
##' @importFrom rlang check_installed
## @importFrom europepmc epmc_hits_trend
##' @importFrom utils modifyList
##' @export
##' @author guangchuang yu
pmcplot <- function(query, period, proportion = TRUE) {
    
    check_installed('europepmc', 'for `pmcplot()`.')
    
    res <- map_df(query, function(x) {
        period <- get("period", parent.env(parent.env(new.env())))	
	y <- europepmc::epmc_hits_trend(query = x, period = period)
        y$query <- x
        return(y)
    })

    mapping <- aes_(x=~year, y = ~query_hits, color = ~query)
    ylab <- "Number of articles"
    if (proportion) {
        mapping <- modifyList(mapping, aes_(y = ~query_hits/all_hits))
        ylab <- "Proportion of articles"
    }
    ggplot(res, mapping) + geom_line() + geom_point() +
        xlab(NULL) + ylab(ylab)
}


