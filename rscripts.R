
library(Hmisc)
setwd("~/phaistos_proj")
dG=read.table("deltaG_new")
deltalogP=read.table("stability_data")

my_frame <- data.frame(dG[,2], deltalogP[,3])
str(my_frame)

hist(deltalogP[,3],probability= TRUE, breaks = 11)


plot(dG[,2], deltalogP[,3],xlab="∆∆G",ylab=expression("∆log(P)"~italic(MuMu)))

rcorr(as.matrix(my_frame), type="pearson")
pname( my_frame)


par(mfrow=c(3,2))
qqnorm(dG[,2])
qqline(dG[,2])

qqnorm(deltalogP[,3])
qqline(deltalogP[,3])

shapiro.test(rnorm(deltalogP[,1]))
shapiro.test(rnorm(dG[,1]))

cor(dG[,2],deltalogP[,3])

library(colorspace)

prop.test(()), correct=FALSE)


barplot(t(as.matrix(dG)), border=NA, main="dG", cex.names=0.6, las=2, ylab="")
barplot(t(as.matrix(deltalogP)), border=NA, main="delta log(P)", cex.names=0.6, las=2, ylab="")