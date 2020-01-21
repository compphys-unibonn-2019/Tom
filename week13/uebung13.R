res <- read.table("corr_boot_sym.tsv")
res_orig <- read.table("orig_corr.tsv")

Tmax <- 48
tmin <- 13
tmax <- 24

fit_data_y <- res_orig[,tmin:tmax]
fit_data_boot_y <- res[,tmin:tmax]
fit_data_x <- tmin:tmax
fit_data_sd <- apply(fit_data_y, 2, sd)

Cov_mat <- cov(fit_data_y)

# Here we try to invert the covariance Matrix. If it is not invertible we perform an uncorrelated fit
inv_Cov <- tryCatch(solve(cov(fit_data_y)), error = function(cond) { 
                    message("Covariance Matrix not invertible!\n")
		    return(diag(1/fit_data_sd^2))
		    }
		    )

x <- fit_data_x
y_mean <- apply(fit_data_y, 2, mean)
#inv_Cov <- diag(1/fit_data_sd^2)### uncomment if you want to fit uncorrelated

model_func1 <- function(x, y, par) {
	Y <- par[1] * exp(-par[2] * x)
	return(Y)
}

model_func2 <- function(x, y, par) {
	Y <- par[1] * exp(-par[2] * x) + par[3] * exp(-par[4] * x)
	return(Y)
}

chisq <- function(x, y, invCov, par) {
	chi <- sum(t(y - par[1] * (exp(-par[2] * x) + exp(-par[2] * (Tmax - x)))) %*% invCov %*% (y - par[1] * (exp(-par[2] * x) + exp(-par[2] * (Tmax - x)))))
	return(chi)
}

###### Fit to the original data
fit_res <- optim(par = c(200, .18), fn = chisq, x = x, y = y_mean, invCov = inv_Cov, method = "BFGS")

dim_data <- dim(fit_data_boot_y)


##### Perform the bootstrap fits
fit_res_boot_vec <- numeric()
for(i in 1:dim_data[1]) {
	y <- as.vector(t(fit_data_boot_y[i,]))
	res_opt <- optim(par = c(200, .18), fn = chisq, x = x, y = y, invCov = inv_Cov, , method = "BFGS")
	fit_res_boot_vec <- c(fit_res_boot_vec, res_opt$par)
}
error_amplitude <- sd(fit_res_boot_vec[seq(1, 2999, 2)])
error_energy <- sd(fit_res_boot_vec[seq(2, 3000, 2)])
fit_res_boot_energy <- array(fit_res_boot_vec, dim = c(2,dim_data[1]))

##### Plotting part
plot(x = x, y = y_mean, xlab = "t", ylab = "C(t)")
x_fit <- seq(tmin, tmax, .1)
y_fit <- fit_res$par[1] * (exp(-fit_res$par[2] * x_fit) + exp(-fit_res$par[2] * (Tmax - x_fit)))

lines(x = x_fit, y = y_fit, col = "red")


##### Result printing part
cat("Fit Result is: c0 = ", fit_res$par[1], "+-", error_amplitude, "\nand E0 =", fit_res$par[2], "+-", error_energy, "\n")
