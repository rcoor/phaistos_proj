library('colorspace')
library('calibrate')
setwd("~/phaistos_proj")

WT=read.table("array_2TRX_wt_mumu.txt", stringsAsFactors=FALSE)
L78K=read.table("array_2TRX_L78K_mumu.txt", stringsAsFactors=FALSE)
L78R=read.table("array_2TRX_L78R_mumu.txt", stringsAsFactors=FALSE)
T77C_V91C=read.table("array_2TRX_T77C_V91C_mumu.txt", stringsAsFactors=FALSE)

par(mfrow=c(1,1))


barplot(t(as.matrix(wt)), names.arg = c(1:77,"*78",79:108), col=rainbow_hcl(2), border=NA, main="MuMu: 2TRX WT", cex.names=0.6, las=2, ylab="log(P)")

barplot(t(as.matrix(L78K)), names.arg = c(1:77,"*78",79:108), col =rainbow_hcl(2), border=NA, main="MuMu: 2TRX L78K", cex.names=0.6, las=2, ylab="log(P)")

barplot(t(as.matrix(L78R)), names.arg = c(1:77,"*78",79:108), col =rainbow_hcl(2), border=NA, main="MuMu: 2TRX L78R", cex.names=0.6, las=2, ylab="log(P)")

barplot(t(as.matrix(T77C_V91C)), names.arg = c(1:76,"C*77",78:90,"C*91",92:108), col =rainbow_hcl(2), border=NA, main="MuMu: 2TRX T77C & V91C", cex.names=0.6, las=2, ylab="log(P)")

log_P_WT <- sum(WT)
log_P_L78K <- sum(L78K)
log_P_L78R <- sum(L78K)
log_P_T77C_V91C <- sum(T77C_V91C)

log_P_WT_L78K <- log_P_WT-log_P_L78K
log_P_WT_L78R <- log_P_WT-log_P_L78R
log_P_WT_T77C_V91C <- log_P_WT-log_P_T77C_V91C

ddG_L78K <- -3.90
ddG_L78R <- -4.0
ddG_T77C_V91C <- 2.1

logP_diff <- c(log_P_WT_L78K,log_P_WT_L78R,log_P_WT_T77C_V91C)
ddG <- c(ddG_L78K,ddG_L78R,ddG_T77C_V91C)
names <- c("L78K","L78R","T77C_V91C")

cov(ddG,logP_diff)
plot(ddG,logP_diff, xlim=c(-5,3), ylim=c(-4,-2.5))
textxy(ddG,logP_diff,names,cex = 0.8)