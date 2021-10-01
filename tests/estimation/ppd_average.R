pacman::p_load("kableExtra", "magrittr", "gtools")
setwd("~/Documents/git/MSc-thesis/code/msMK")
source("tests/tests.R")
devtools::load_all()

setwd("/home/dede/Documents/git/MSc-thesis/code/tests/estimation")
# 2 dimensions
dirs = list.dirs(recursive = FALSE)
dirs = dirs[grepl("_2", dirs)]
for(w in 1:length(dirs)){
    imgname = paste0(dirs[w], "/img/", "average_ppdSmooth-indTRUE-red.png")
    if(!file.exists(imgname)){
        pred_list = vector(mode = "list", length = 10)
        cat(dirs[w], "\n")
        nsim = 50
        for(j in 1:nsim){
            cat("j:", j, "\n")
            load(paste0(dirs[w], "/i", j,"-indTRUE.Rdata"))
            pred_list[[j]] = predictive_smooth_average(mcmc.msmk, y.stand, ylim = c(-4.5,4.5))
        }
        save(pred_list, file = paste0(dirs[w], "pred_list-indTRUE.Rdata"))
        pred = matrix(0, nrow = NROW(pred_list[[1]]$fhat), ncol = NCOL(pred_list[[1]]$fhat))
        for(j in 1:nsim){
            pred = pred + pred_list[[j]]$fhat
        }
        pred = pred / nsim

        png(imgname, width = 2200, height = 1200, res = 100)
        x = pred_list[[1]]$x1
        y = pred_list[[1]]$x2
        contour(x, y, pred)
        densContour(df, xlim = c(min(x), max(x)), ylim = c(min(y), max(y)), nlevels=20, col = "darkorange",
        add = TRUE)
        dev.off()
    }
}


for(w in 1:length(dirs)){
    imgname = paste0(dirs[w], "/img/", "average_ppdSmooth-indFALSE-red.png")
    if(!file.exists(imgname)){
        pred_list = vector(mode = "list", length = 10)
        cat(dirs[w], "\n")
        nsim = 50
        for(j in 1:nsim){
            cat("j:", j, "\n")
            load(paste0(dirs[w], "/i", j,"-indFALSE.Rdata"))
            pred_list[[j]] = predictive_smooth_average(mcmc.msmk, y.stand, ylim = c(-4.5,4.5))
        }
        save(pred_list, file = paste0(dirs[w], "pred_list-indFALSE.Rdata"))
        pred = matrix(0, nrow = NROW(pred_list[[1]]$fhat), ncol = NCOL(pred_list[[1]]$fhat))
        for(j in 1:nsim){
            pred = pred + pred_list[[j]]$fhat
        }
        pred = pred / nsim
        
        png(imgname, width = 2200, height = 1200, res = 100)
        x = pred_list[[1]]$x1
        y = pred_list[[1]]$x2
        contour(x, y, pred)
        densContour(df, xlim = c(min(x), max(x)), ylim = c(min(y), max(y)), nlevels=20, col = "darkorange",
        add = TRUE)
        dev.off()
    }
}


# for(w in 1:length(dirs)){
#     imgname = paste0(dirs[w], "/img/", "expectedScale-indTRUE.png")
#     load(paste0(dirs[w], "/i", 5,"-indTRUE.Rdata"))
#     png(imgname, width = 2200, height = 1200, res = 100) 
#     msMK.plot.depth(mcmc.msmk, TRUE)
#     dev.off()
#     imgname = paste0(dirs[w], "/img/", "expectedScale-indFALSE.png")
#     load(paste0(dirs[w], "/i", 5,"-indFALSE.Rdata"))
#     png(imgname, width = 2200, height = 1200, res = 100) 
#     msMK.plot.depth(mcmc.msmk, TRUE)
#     dev.off()
# }
