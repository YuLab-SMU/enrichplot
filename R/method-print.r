##' @method print enrichplotDot
##' @export
print.enrichplotDot <- function(x, ...) {
    p <- ggfun::set_point_legend_shape(x)
    class(p) <- class(p)[-1]
    print(p)
}
