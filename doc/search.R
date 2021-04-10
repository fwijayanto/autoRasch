## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(autoRasch)

## ------------------------------------------------------------------------
ipoqll_score <- compute_score(short_poly_data, incl_set = c(7:9), type = "ipoqll")
summary(ipoqll_score)

## ------------------------------------------------------------------------
ipoqll_scores <- compute_scores(short_poly_data, incl_sets = rbind(c(1:3),c(7:9)), type = "ipoqll", cores = 2)
View(ipoqll_scores[,1:12])

## ------------------------------------------------------------------------
stepwise_res <- stepwise_search(short_poly_data, criterion = "ipoqll", cores = 2, 
                                isTracked = TRUE)

## ----fig.height=3.5, fig.width=7-----------------------------------------
plot_search(stepwise_res, type="l")

