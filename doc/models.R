## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(autoRasch)

## ------------------------------------------------------------------------
pcm_res <- pcm(polydif_inh_dset[,13:19])

## ------------------------------------------------------------------------
pcmdif_res <- pcm_dif(polydif_inh_dset[,13:19], groups_map = c(rep(0,245),rep(1,245)))

## ------------------------------------------------------------------------
summary(pcm_res, par="beta")

## ------------------------------------------------------------------------
summary(pcmdif_res, par="delta")

## ----fig.height=3.5, fig.width=7-----------------------------------------
plot_PImap(pcm_res, main = "Person-Item map of the PCM")

## ----fig.height=3.5, fig.width=7-----------------------------------------
plot_PImap(pcmdif_res, main = "Person-Item map of the PCM-DIF")

## ----fig.height=3.5, fig.width=7-----------------------------------------
plot_ICC(pcm_res, itemno = 5, main = "ICC of I 17; estimated using PCM")

## ------------------------------------------------------------------------
pcm_fit <- fitStats(pcm_res)
itemfit(pcm_fit)

## ------------------------------------------------------------------------
pcmdif_fit <- fitStats(pcmdif_res)
itemfit(pcmdif_fit)

## ------------------------------------------------------------------------
gpcm_res <- gpcm(polydif_inh_dset[,13:19])

## ------------------------------------------------------------------------
gpcmdif_res <- gpcm_dif(polydif_inh_dset[,13:19], groups_map = c(rep(0,245),rep(1,245)))

## ------------------------------------------------------------------------
summary(pcm_res, par="alpha")

## ------------------------------------------------------------------------
summary(pcmdif_res, par="delta")

