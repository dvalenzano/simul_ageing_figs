dev.off()

pop <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/popres2.csv", header=T)
setwd('/Volumes/group_dv/personal/DValenzano/papers/simulation_arXiv/Figure3/')

dev.off()
pdf("Figure3A.pdf", width=4, height=3.2)
fig3 <- read.csv('sc_25k.csv', sep=',', header=F)
x <- c(1:5000)
y1 <- fig3[1:(length(fig3)-20001)]
y2 <- rep(5000, 5000)
plot(x,y1, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,20000), pin=c(3,2))
lines(x, y2, col="red", lwd=2)
dev.off()


pdf("Figure3B.pdf", width=4, height=3.2)
x <- pop$time
y2 <- pop$resources
y1 <- pop$population
plot(x,y2, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,20000), pin=c(3,2))
lines(x, y1, col="blue", lwd=2)
dev.off()

