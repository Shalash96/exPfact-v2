library(mclust)

args <- commandArgs(trailingOnly = TRUE)
XXX <- as.integer(args[1])
YYY <- as.integer(args[2])
sp_path <- if (length(args) >= 3) args[3] else "all.sp"

data <- read.table(sp_path)
out <- paste0(XXX, "-", YYY)

fit <- Mclust(data[, XXX:YYY], G = 1:99, modelNames = "VVV")
write.table(file = paste0(out, ".mod"), fit$parameters$mean)
write.table(file = paste0(out, ".pro"), fit$parameters$pro)

dat <- numeric()
for (i in 1:length(fit$parameters$pro)) {
  dat <- cbind(dat, sqrt(diag(fit$parameters$variance$sigma[,,i])))
}
write.table(file = paste0(out, ".var"), dat)
write.table(file = paste0(out, ".bic"), fit$BIC)
write.table(file = paste0(out, ".z"), fit$z)
