library(plotrix)
lp <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/large-sex-popres.csv", header=T)
lc <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/large-sex-constres.csv", header=T)
sp <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-sex-popres.csv", header=T)
sc <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-sex-constres.csv", header=T)

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6A.pdf", width=4, height=3.2)
mfrow=c(1,1)
lpe = lp[1:5001,]
plot(lpe$time,lpe$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lpe$time, lpe$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6B.pdf", width=4, height=3.2)
mfrow=c(1,1)
lpl = lp[7001:12000,]
plot(lpl$time,lpl$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lpl$time, lpl$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6C.pdf", width=4, height=3.2)
mfrow=c(1,1)
lce = lc[1:5001,]
plot(lce$time,lce$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lce$time, lce$resources, col="red", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6D.pdf", width=4, height=3.2)
mfrow=c(1,1)
lcl = lc[7001:12000,]
plot(lcl$time,lcl$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lcl$time, lcl$resources, col="red", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6E.pdf", width=4, height=3.2)
mfrow=c(1,1)
spe = sp[1:5001,]
plot(spe$time,spe$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,6000), pin=c(3,2))
lines(spe$time, spe$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6F.pdf", width=4, height=3.2)
mfrow=c(1,1)
spl = sp[7001:12000,]
plot(spl$time,spl$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,6000), pin=c(3,2))
lines(spl$time, spl$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6G.pdf", width=4, height=3.2)
mfrow=c(1,1)
sce = sc[1:5001,]
plot(sce$time,sce$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,6000), pin=c(3,2))
lines(sce$time, sce$resources, col="red", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Figure6H.pdf", width=4, height=3.2)
mfrow=c(1,1)
scl = sc[7001:12000,]
plot(scl$time,scl$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,6000), pin=c(3,2))
lines(scl$time, scl$resources, col="red", lwd=2)
dev.off()

