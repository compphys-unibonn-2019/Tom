library(hadron)

filename <- "/home/matthias/Masterarbeit/Daten/pi_corr_p0_A8024.dat" 
corr <- readtextcf(filename, ind.vector = c(2, 2), T=48) #read file with corr fct data
corr$icf <- corr$icf * 0 #set imaginary part of C(t) to 0
corr_boot <- bootstrap.cf(corr, boot.R = 1500, boot.l = 1, seed = 1)
em <- bootstrap.effectivemass(corr_boot, type = "acosh")

plot(corr_boot, main = "Pion Correlator A80.24", xlab = "t", ylab = "C(t)")
dev.new()
plot(em, main = "Effective Mass for a Pion in A80.24 at P^2=0", xlab = "t", ylab = expression(m[eff](t)))
