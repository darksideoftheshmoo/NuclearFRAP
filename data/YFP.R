if(is.element("d",ls())) warning("d was replaced",call. = FALSE)

load("YFPdata.rda")
suppressWarnings(d<-update_img.path(d,system.file('img', package='NuclearFRAP')))
