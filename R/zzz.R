globals <- new.env()

.onLoad <- function(libname, pkgname) {
  globals$glab_brand <- c("AZ", "PF")
  globals$glab_trt <- rbind(c("CVD+FLU", "CVD+PBO"),
                            c("PBO", "FLU"))
}

