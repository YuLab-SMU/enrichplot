##' ggplot theme of DOSE
##'
##' @title theme_dose
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @export
##' @param font.size font size
theme_dose <- function(font.size=14) {
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black",
                                     size = font.size, vjust =1 ),
          axis.text.y = element_text(colour = "black",
                                     size = font.size, hjust =1 ),
          axis.title = element_text(margin=margin(10, 5, 0, 0),
                                    color = "black",
                                    size = font.size),
          axis.title.y = element_text(angle=90)
          )
}
