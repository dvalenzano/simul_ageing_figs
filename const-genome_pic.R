##GOAL: To generate the simulation_arXiv figures relative to S, R throughout the simulation, standard deviation and s_i, r_i

#### First, one panel at the time ####

abcdz <- read.csv('/Volumes/group_dv/personal/DValenzano/papers/simulation_arXiv/Figure3/first-run/first2runs.csv', sep=',', head=T)
xabcdz <- c(1:(length(abcdz$het)/5))
abcdz0 <- subset(abcdz, group == 0)
abcdz1 <- subset(abcdz, group == 1)
abcdz5k <- subset(abcdz, group == 3)
abcdz10k <- subset(abcdz, group == 6)
abcdz60k <- subset(abcdz, group == "z")

x2 <- c(1:71, 16:70)

## AS A LINE PLOT ##
pdf("Figure4_const.pdf", width=4, height=3.2)
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
sd_abcdz <- read.csv('/Volumes/group_dv/personal/DValenzano/papers/simulation_arXiv/Figure3/het-sd.csv', sep=',', head=T)
sd_xabcdz <- c(1:(length(sd_abcdz$het)/5))
sd_abcdz0 <- subset(sd_abcdz, group == 0)
sd_abcdz1 <- subset(sd_abcdz, group == 1)
sd_abcdz5k <- subset(sd_abcdz, group == 3)
sd_abcdz10k <- subset(sd_abcdz, group == 6)
sd_abcdz60k <- subset(sd_abcdz, group == "z")


pdf(file="sd-het_const.pdf", width=4, height=3.2)
plot(sd_xabcdz, sd_abcdz60k$het, ylim = c(0.75, 1.2), ylab = "Genetic Variance", xlab = "Age", type="l", lwd=1.5, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

### THIS IS ADDED ON 16-SEP-2015 #######
setwd('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/')

abcd <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-large_surv.csv', sep=',', head=T)
xabcd <- c(1:(length(abcd$group)/4))
abcd_lps <- subset(abcd, group == 'lps')
abcd_lcs <- subset(abcd, group == 'lcs')
abcd_sps <- subset(abcd, group == 'sps')
abcd_scs <- subset(abcd, group == 'scs')
x2 <- c(1:71, 16:70)

## AS A LINE PLOT ##
mfrow=c(1,1)
pdf("Figure6B_pop.pdf", width=4, height=3.2)
plot(xabcd, abcd_lps$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, abcd_sps$surv, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

mfrow=c(1,1)
pdf("Figure6C_const.pdf", width=4, height=3.2)
plot(xabcd, abcd_lcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, abcd_scs$surv, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

### THIS IS ADDED ON 17-SEP-2015 ### I AM ADDING THE 'JUNK' DNA LEVEL ####

Jnk <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Jnk_sex.csv', sep=',', head=T)
Jnk <- Jnk[1:126,]

## AS A LINE PLOT ##
mfrow=c(1,1)
pdf("Figure6B2_pop.pdf", width=4, height=3.2)
plot(xabcd, abcd_lps$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, abcd_sps$surv, lty=1, lwd=3, col="red")
lines(xabcd, Jnk$Jsp, lty=2, lwd=3, col="red")
lines(xabcd, Jnk$Jlp, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

mfrow=c(1,1)
pdf("Figure6C2_const.pdf", width=4, height=3.2)
plot(xabcd, abcd_lcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, abcd_scs$surv, lty=1, lwd=3, col="red")
lines(xabcd, Jnk$Jsc, lty=2, lwd=3, col="red")
lines(xabcd, Jnk$Jlc, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

######NOW FOR FIGURE 4#########

pdf("Figure4_new_const.pdf", width=4, height=3.2)
plot(xabcdz, abcdz0$het, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcdz, abcdz1$het, col=2, lwd=3)
lines(xabcdz, abcdz5k$het, col=3, lwd =3)
lines(xabcdz, abcdz10k$het, col=4, lwd=3)
lines(xabcdz, abcdz60k$het, col=5, lwd=3)
lines(xabcd, Jnk$Jlc, lty=2, lwd=3, col=5)
abline(v=16, lwd=2, lty=2)
abline(v=72, lwd=4, lty=1)
dev.off()

########### ASEX - SMALL, LARGE, VERY LARGE POP #############

bde_as <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/const-small-large_surv_repr.csv', sep=',', head=T)
xbde <- c(1:(length(bde_as$group)/3))
bde_lcs <- subset(bde_as, group == 'lcs')
bde_scs <- subset(bde_as, group == 'scs')
bde_llcs <- subset(bde_as, group == 'llcs')
x2 <- c(1:71, 16:70)

## AS A LINE PLOT ##
mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_pop.pdf", width=4, height=3.2)
plot(xbde, bde_lcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, col="gray48", lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xbde, bde_scs$surv, lty=1, lwd=3, col="red")
lines(xbde, bde_llcs$surv, lty=1, lwd=3, col="blue")
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

# Now I add the junk DNA 

Jnk_as <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Jnk_sex-repr_as.csv', sep=',', head=T)
Jnk_as <- Jnk_as[1:126,]

## AS A LINE PLOT ##
mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_lc.pdf", width=4, height=3.2)
plot(xbde, bde_lcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, Jnk_as$Jlc, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_sc.pdf", width=4, height=3.2)
plot(xbde, bde_scs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, Jnk_as$Jsc, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_llc.pdf", width=4, height=3.2)
plot(xbde, bde_llcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcd, Jnk_as$Jllc, lty=2, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
dev.off()

########## 20 SEP 2015 ############
### ASEX POP LARGE AND SMALL ######
###################################

abcde_as <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/const-small-large_surv_repr.csv', sep=',', head=T)
xabcde <- c(1:(length(abcde_as$group)/5))
abcde_lcs <- subset(abcde_as, group == 'lcs')
abcde_lps <- subset(abcde_as, group == 'lps')
abcde_scs <- subset(abcde_as, group == 'scs')
abcde_sps <- subset(abcde_as, group == 'sps')
abcde_llcs <- subset(abcde_as, group == 'llcs')
x2 <- c(1:71, 16:70)

Jnk_as <- read.csv('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Jnk_sex-repr_as.csv', sep=',', head=T)
Jnk_as <- Jnk_as[1:126,]

y.high_c<- rep(max(Jnk_as$Jlc, Jnk_as$Jsc, Jnk_as$Jllc), length(xabcde))
y.low_c<- rep(min(Jnk_as$Jlc, Jnk_as$Jsc, Jnk_as$Jllc), length(xabcde))

library(grDevices)

mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_slllc.pdf", width=4, height=3.2)
plot(xabcde, abcde_lcs$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcde, abcde_scs$surv,lty=1, lwd=3, col="red")
lines(xabcde, abcde_llcs$surv,lty=2, lwd=3)
#lines(xabcde, Jnk_as$Jsc, lty=1, lwd=3, col="red")
#lines(xabcde, Jnk_as$Jlc, lty=1, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
polygon(c(xabcde,rev(xabcde)), c(y.high_c, rev(y.low_c)), border=NA, col=rgb(0.41,0.41,0.41,0.5) )
dev.off()

y.high_p<- rep(max(Jnk_as$Jlp, Jnk_as$Jsp), length(xabcde))
y.low_p<- rep(min(Jnk_as$Jlp, Jnk_as$Jsp), length(xabcde))

mfrow=c(1,1)
pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7D_slp.pdf", width=4, height=3.2)
plot(xabcde, abcde_lps$surv, ylim = c(0, 1), ylab = "Frequency of 1s / S/R_i", xlab = "Age", type="l", lwd=3, lty=1, xaxt="n", bty="n")
axis(1, at=1:126, labels=x2, tick=T, lwd.ticks=0)
lines(xabcde, abcde_sps$surv,lty=1, lwd=3, col="red")
#lines(xabcde, Jnk_as$Jsc, lty=1, lwd=3, col="red")
#lines(xabcde, Jnk_as$Jlc, lty=1, lwd=3)
abline(v=16, lwd=1, lty=2)
abline(v=72, lwd=2, lty=1)
polygon(c(xabcde,rev(xabcde)), c(y.high_p, rev(y.low_p)), border=NA, col=rgb(0.41,0.41,0.41,0.5) )
dev.off()


