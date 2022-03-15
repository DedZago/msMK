pacman::p_load("kableExtra", "magrittr", "gtools", "stringr")
setwd("/home/dede/Documents/git/msMK/")
source("tests/tests.R")
devtools::load_all()

setwd("/home/dede/Documents/git/msMK/tests/estimation")
# 2 dimensions
dirs = list.dirs(recursive = FALSE)
dirs = dirs[grepl("_2", dirs)]
for(w in 1:length(dirs)){
    imgname = paste0(dirs[w], "/img/", str_replace(dirs[w], "./r", "f"), ".png")
    if(!file.exists(imgname)){
        dir.create(paste0(dirs[w], "/img"), showWarnings = FALSE)
        df = get(str_replace(dirs[w], "./r", "f"))
        cat(dirs[w], "\n")
        png(imgname)
        densContour(df, xlim = c(-4,4), ylim = c(-4,4), nlevels=160, col = "black")
        dev.off()
    }
}
