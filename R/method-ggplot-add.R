##' @importFrom ggplot2 ggplot_add
##' @method ggplot_add autofacet
##' @export
ggplot_add.autofacet <- function(object, plot, object_name) {
    d <- plot$data
    nn <- names(d)
    if ('category' %in% nn) {
        var <- "category"
    } else if ('ONTOLOGY' %in% nn) {
        var <- 'ONTOLOGY'
    } else {
        message("not supported")
        return(plot)
    }

    if (!is.null(object$levels)) {
        d[[var]] <- factor(d[[var]], levels = object$levels)
        plot$data <- d
    }
    if (object$by == 'row') {
        obj <- facet_grid(.data[[var]] ~ ., scales=object$scales)
    } else {
        obj <- facet_grid(. ~ .data[[var]], scales=object$scales)
    }
    ggplot_add(obj, plot, object_name)
}
