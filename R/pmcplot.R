#' PubMed Central Trend plot
#'
#'
#' @title pmcplot
#' @param query query terms
#' @param period period of query in the unit of year
#' @param proportion If TRUE, use query_hits/all_hits, otherwise use query_hits.
#' @param data Optional data.frame with `query`, `year`, `query_hits`, and
#'   `all_hits` columns. When supplied, no Europe PMC request is made.
#' @return ggplot object
#' @importFrom purrr map_df
#' @importFrom rlang check_installed
#' @importFrom utils modifyList
#' @export
#' @author Guangchuang Yu
pmcplot <- function(query, period, proportion = TRUE, data = NULL) {
    if (is.null(data)) {
        require_suggested('europepmc', 'for `pmcplot()`.')

        res <- map_df(query, function(x) {
            period <- get("period", parent.env(parent.env(new.env())))
            y <- europepmc::epmc_hits_trend(query = x, period = period)
            y$query <- x
            return(y)
        })
    } else {
        required <- c("query", "year", "query_hits", "all_hits")
        if (!is.data.frame(data) || !all(required %in% colnames(data))) {
            stop(
                "`data` must contain columns: ",
                paste(required, collapse = ", ")
            )
        }
        res <- data[
            data$query %in% query & data$year %in% period,
            required,
            drop = FALSE
        ]
        if (!nrow(res)) {
            stop("`data` has no rows matching `query` and `period`.")
        }
    }

    res$query <- factor(res$query, levels = query)
    mapping <- aes(x = .data$year, y = .data$query_hits, color = .data$query)
    ylab <- "Number of articles"
    if (proportion) {
        mapping <- modifyList(mapping, aes(y = .data$query_hits / .data$all_hits))
        ylab <- "Proportion of articles"
    }
    ggplot(res, mapping) +
        geom_line(linewidth = 0.6) +
        geom_point(size = 1.5) +
        scale_x_continuous(breaks = period) +
        xlab(NULL) +
        ylab(ylab)
}


