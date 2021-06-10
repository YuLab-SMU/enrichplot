##' Get the frequency of each word in a vector of terms.
##' 
##' @param wordd clusters, a vector of terms.
##' @noRd
get_word_freq <- function(wordd){     
    dada <- strsplit(wordd, " ")
    didi <- table(unlist(dada))
    didi <- didi[order(didi, decreasing = TRUE)]
    # Get the number of each word
    word_name <- names(didi)
    fun_num_w <- function(ww){
        sum(vapply(dada, function(w){ww %in% w}, FUN.VALUE = 1))
    }
    word_num <- vapply(word_name, fun_num_w, FUN.VALUE = 1)
    word_w <- word_num[order(word_num, decreasing = TRUE)]
}

##' Use wordcloud algorithm to get group tags
##' 
##' @param cluster a cluster name
##' @param ggData the data section of the ggraph object, 
##' which contains clustering information.
##' @param nWords the number of words in the cluster tags 
##' @importFrom magrittr %>%
##' @noRd
get_wordcloud <- function(cluster, ggData, nWords){
    words <- ggData$name %>%
        gsub(" in ", " ", .) %>%
        gsub(" [0-9]+ ", " ", .) %>%
        gsub("^[0-9]+ ", "", .) %>%
        gsub(" [0-9]+$", "", .) %>%
        gsub(" [A-Za-z] ", " ", .) %>%
        gsub("^[A-Za-z] ", "", .) %>%
        gsub(" [A-Za-z]$", "", .) %>%
        gsub(" / ", " ", .) %>%
        gsub(" and ", " ", .) %>%
        gsub(" of ", " ", .) %>%
        gsub(",", " ", .) %>%
        gsub(" - ", " ", .)
    net_tot <- length(words)

    clusters <- unique(ggData$color2)
    words_i <- words[which(ggData$color2 == cluster)]

    sel_tot <- length(words_i)
    sel_w <- get_word_freq(words_i)
    net_w_all <- get_word_freq(words)
    net_w <- net_w_all[names(sel_w)]
    tag_size <- (sel_w/sel_tot)/(net_w/net_tot)
    tag_size <- tag_size[order(tag_size, decreasing = TRUE)]
    nWords <- min(nWords, length(tag_size))
    tag <- names(tag_size[seq_len(nWords)])

    # Order of words
    dada <- strsplit(words_i, " ")
    len <- vapply(dada, length, FUN.VALUE=1)
    rank <- NULL
    for(i in seq_len(sel_tot)) {
        rank <- c(rank, seq_len(len[i]))
    }

    word_data <- data.frame(word = unlist(dada), rank = rank)
    word_rank1 <- stats::aggregate(rank ~ word, data = word_data, sum)
    rownames(word_rank1) <- word_rank1[, 1]

    word_rank1 <- word_rank1[names(sel_w), ]
    # Get an average ranking order
    word_rank1[, 2] <- word_rank1[, 2]/as.numeric(sel_w)
    tag_order <- word_rank1[tag, ]
    tag_order <- tag_order[order(tag_order[, 2]), ]
    tag_clu_i <- paste(tag_order$word, collapse=" ")
}













