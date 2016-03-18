# This script plots OU-BDSKY plot

#  origin_post is a posterior vector of origins
#  r0 is a data table of posterior samples of r0 vectors, one row per sample
#  time_grid is a vector of times to evaluate the skyline at
bdsky_post <- function(origin_post, r0, time_grid) {
	
	r0_time_gridded <- list()
	
	n <- ncol(r0)
		
	for (s in 1:length(origin_post)) {
		
		origin <- origin_post[s]
		r0_vec <- r0[s,]
				
		ind <- pmax(1,n - floor(time_grid / origin * n))
		r0_time_gridded[[s]] <- r0_vec[ind]			
	}
	
	return (r0_time_gridded)
}

setwd("~/Git/bdsky/examples/")

lf <- read.table("HCV_oup_40_1447852063188.log", sep="\t", header=T)

origin_post <- lf$orig_root

r0_subset <- lf[grepl("R0", names(lf))]

time_grid <- 1:400

bdskypost <- bdsky_post(origin_post, r0_subset, time_grid)

pdf("hcv_oubdsky.pdf")
plot(time_grid,bdskypost[[950]],type='S', xlab="Time (years before present)", ylab="R0",col=rgb(0,0,1,0.1))
for (s in 20:200*50) {
	lines(time_grid, bdskypost[[s]], type='S',col=rgb(0,0,1,0.1))
}
dev.off()
