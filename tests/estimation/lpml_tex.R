pacman::p_load("kableExtra", "magrittr", "gtools")

# 2 dimensions
dirs = list.dirs(recursive = FALSE)
dirs = dirs[grepl("_2", dirs)]
dirs = mixedsort(dirs)  
lpml_list = vector(mode = "list", length = length(dirs)*3)
setwd("/home/dede/Documents/git/msMK/tests/estimation")

ns = c(250, 500, 1000)
for(n in 1:length(ns)){
    for(i in 1:length(dirs)){
        print((n-1)*length(dirs) + i)
        load(paste0(dirs[i], "/lpml-indTRUE", "-n", ns[n], ".Rdata"))
        temp = lpml_mat
        load(paste0(dirs[i], "/lpml-indFALSE", "-n", ns[n], ".Rdata"))
        temp = cbind(temp, lpml_mat[ , 1 ])
        load(paste0(dirs[i], "/lpml-PYfix", "-n", ns[n],".Rdata"))
        temp = cbind(temp, lpml_mat[ , 1 ])
        temp = temp[ , c(1, 3, 4, 2) ]
        colnames(temp) = c("Product", "Full covariance", "PY", "True")
        # temp2 = temp[ , c(1,2,3) ] / temp[, 4]
        # temp2 = (temp2)^-1
        # print(temp2)
        # colnames(temp2) = c("Product", "Full covariance", "PY")
        lpml_mean = apply(na.omit(temp), 2, mean) %>% round(3)
        lpml_sd = apply(na.omit(temp), 2, sd) %>% round(3)
        lpml_list[[(n-1)*length(dirs) + i]] = paste0(format(lpml_mean, nsmall = 1), " (", format(lpml_sd, nsmall=1), ")")
    }
}
do.call(rbind, lpml_list) %>% 
    set_rownames(rep(paste0("(", letters[1:length(dirs)], ")"), 3))%>% 
    set_colnames(c("Product", "Full covariance", "PY", "True")) %>% 
    kbl(format="latex", booktabs = TRUE, linesep = "") %>% 
    group_rows(index = c("$n = 250$" = 6, "$n = 500$" = 6, "$n = 1000$" = 6), escape = FALSE, bold = FALSE)

# TODO: Bold for better model?

# Higher dimensions

dirs = list.dirs(recursive = FALSE)
notIdx = grepl("_2", dirs)
dirs = dirs[!notIdx]
first = grepl("_5", dirs)
d = sum(first)
last = grepl("_10", dirs)
beg = mixedsort(dirs[first])
end = mixedsort(dirs[last])
dirs = c(beg, end)

lpml_list = vector(mode = "list", length = length(dirs))

for(i in 1:length(dirs)){
    cat(dirs[i], "TRUE\n")
    load(paste0(dirs[i], "/lpml-indTRUE-n1000.Rdata"))
    temp = lpml_mat
    cat(dirs[i], "FALSE\n")
    load(paste0(dirs[i], "/lpml-indFALSE-n1000.Rdata"))
    temp = cbind(temp, lpml_mat[ , 1 ])
    load(paste0(dirs[i], "/lpml-PYfix-n1000.Rdata"))
    temp = cbind(temp, lpml_mat[ , 1 ])
    temp = temp[ , c(1, 3, 4, 2) ]
    colnames(temp) = c("Product", "Full covariance", "DP", "True")
    temp2 = temp[ , c(1,2,3) ] / temp[ , 4 ]
    temp2 = (temp2)^-1
    colnames(temp2) = c("Product", "Full covariance", "DP")
    lpml_mean = apply(na.omit(temp2), 2, mean) %>% round(3)
    lpml_sd = apply(na.omit(temp2), 2, sd) %>% round(3)
    
    lpml_list[[i]] = paste0(format(lpml_mean, nsmall = 1), " (", format(lpml_sd, nsmall=1), ")")
}
do.call(rbind, lpml_list) %>% 
    set_rownames(rep(paste0("(", letters[1:d], ")"), length(dirs)/d))%>% 
    set_colnames(c("Product", "Full covariance", "PY")) %>% 
    kbl(format="latex", booktabs = TRUE) %>% 
    group_rows(index = c("$d = 5$" = 6, "$d = 10$" = 6), escape = FALSE, bold = FALSE)
