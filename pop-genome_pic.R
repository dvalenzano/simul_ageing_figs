##GOAL: To generate the simulation_arXiv figures relative to S, R throughout the simulation, standard deviation and s_i, r_i
setwd('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/')

#### First, one panel at the time ####

abcdz <- read.csv('./het-freq.csv', sep=',', head=T)
xabcdz <- c(1:(length(abcdz$het)/5))
abcdz0 <- subset(abcdz, group == 0)
abcdz1 <- subset(abcdz, group == 1)
abcdz5k <- subset(abcdz, group == 3)
abcdz10k <- subset(abcdz, group == 6)
abcdz60k <- subset(abcdz, group == "z")

x2 <- c(1:71, 16:70)

## AS A LINE PLOT ##
pdf("Figure4_pop.pdf", width=4, height=3.2)
plot(xabcdz, abcdz0$het, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcdz, abcdz1$het, col=2, lwd=3)
lines(xabcdz, abcdz5k$het, col=3, lwd =3)
lines(xabcdz, abcdz10k$het, col=4, lwd=3)
lines(xabcdz, abcdz60k$het, col=5, lwd=3)
abline(v=16, lwd=2, lty=2)
abline(v=72, lwd=4, lty=1)
dev.off()

###### PLOTTING STANDARD DEVIATION ######
sd_abcdz <- read.csv('./het-sd.csv', sep=',', head=T)
sd_xabcdz <- c(1:(length(sd_abcdz$het)/5))
sd_abcdz0 <- subset(sd_abcdz, group == 0)
sd_abcdz1 <- subset(sd_abcdz, group == 1)
sd_abcdz5k <- subset(sd_abcdz, group == 3)
sd_abcdz10k <- subset(sd_abcdz, group == 6)
sd_abcdz60k <- subset(sd_abcdz, group == "z")

pdf(file="sd-het_pop.pdf", width=4, height=3.2)
plot(sd_xabcdz, sd_abcdz60k$het, ylim = c(0.75, 1.2), ylab = "Genetic Variance", xlab = "Age", type="l", lwd=1.5, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

### THIS IS ADDED ON 17-SEP-2015 ### I AM ADDING THE 'JUNK' DNA LEVEL ####

Jnk <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Jnk_sex', sep=',', head=T)
Jnk <- Jnk[1:126,]

pdf("Figure4_pop_new.pdf", width=4, height=3.2)
plot(xabcdz, abcdz0$het, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcdz, abcdz1$het, col=2, lwd=3)
lines(xabcdz, abcdz5k$het, col=3, lwd =3)
lines(xabcdz, abcdz10k$het, col=4, lwd=3)
lines(xabcdz, abcdz60k$het, col=5, lwd=3)
lines(xabcd, Jnk$Jlp, lty=2, lwd=3, col=5)
abline(v=16, lwd=2, lty=2)
abline(v=72, lwd=4, lty=1)
dev.off()


