lp_as <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/large-sex-popres.csv", header=T)
lc_as <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/large-sex-constres.csv", header=T)

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7A1.pdf", width=4, height=3.2)
mfrow=c(1,1)
lpe_as = lp_as[1:5001,]
plot(lpe_as$time,lpe_as$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lpe_as$time, lpe_as$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7A2.pdf", width=4, height=3.2)
mfrow=c(1,1)
lpl_as = lp_as[7001:12000,]
plot(lpl_as$time,lpl_as$resources, type = 'l', col="red", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lpl_as$time, lpl_as$population, col="blue", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7A3.pdf", width=4, height=3.2)
mfrow=c(1,1)
lce_as = lc_as[1:5001,]
plot(lce_as$time,lce_as$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lce_as$time, lce_as$resources, col="red", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7A4.pdf", width=4, height=3.2)
mfrow=c(1,1)
lcl_as = lc_as[7001:12000,]
plot(lcl_as$time,lcl_as$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,22000), pin=c(3,2))
lines(lcl_as$time, lcl_as$resources, col="red", lwd=2)
dev.off()

##VERY-LARGE POPULATION, CONSTANT CONDITIONS###

llc_as <- read.csv("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/very-large-sex-constres.csv", header=T)

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7C1.pdf", width=4, height=3.2)
mfrow=c(1,1)
llce_as = llc_as[1:5001,]
plot(llce_as$time,llce_as$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,38000), pin=c(3,2))
lines(llce_as$time, llce_as$resources, col="red", lwd=2)
dev.off()

pdf("/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Figure7C2.pdf", width=4, height=3.2)
mfrow=c(1,1)
llcl_as = llc_as[7001:12000,]
plot(llcl_as$time,llcl_as$population, type = 'l', col="blue", lwd = 1, xlab="Stage", ylab="N", bty="n", ylim=c(0,38000), pin=c(3,2))
lines(llcl_as$time, llcl_as$resources, col="red", lwd=2)
dev.off()

#############################
