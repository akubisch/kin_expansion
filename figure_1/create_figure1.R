
shuf<-read.table("shuffle/data/output_shuffle.txt",header=T)
con<-read.table("control/data/output_noshuffle.txt",header=T)

shuf_m_fst <- tapply(shuf$fst_mar,list(shuf$t),mean)
con_m_fst <- tapply(con$fst_mar,list(con$t),mean)

shuf_m_r <- tapply(shuf$real_b,list(shuf$t),median)
con_m_r <- tapply(con$real_b,list(con$t),median)

shuf_m_e <- tapply(shuf$emi_mar,list(shuf$t),median)
con_m_e <- tapply(con$emi_mar,list(con$t),median)

postscript("figure1.eps",width=3,height=9,paper="special",horizontal=F,title="Kubisch et al. - Figure 1",pointsize=12)

xr<-c(1000,1500)

par(mfrow=c(3,1),mar=c(3,5,1,1),oma=c(2,0,0,0))

plot(shuf_m_r~unique(shuf$t),pch=16,cex=.95,ylim=c(50,200),#ylim=range(c(shuf2$emi_mar,noshuf2$emi_mar)),
     ylab="range border position",bty="l",cex.axis=1.5,cex.lab=1.85,xlim=xr,xlab="",type="l",lwd=2,col="grey50",lty="dashed")
lines(con_m_r~unique(con$t),lwd=2)
text(1030,par("usr")[[4]]*0.95,"A",cex=2)

plot(shuf_m_e~unique(shuf$t),pch=16,cex=.95,ylim=c(0,0.45),#ylim=range(c(shuf2$emi_mar,noshuf2$emi_mar)),
     ylab="emigration rate",bty="l",cex.axis=1.5,cex.lab=1.85,xlim=xr,xlab="",type="l",lwd=2,col="grey50",lty="dashed")
lines(con_m_e~unique(con$t),lwd=2)
text(1030,par("usr")[[4]]*0.95,"B",cex=2)

plot(shuf_m_fst~unique(shuf$t),pch=16,cex=.95,ylim=c(0,0.65),#ylim=range(c(shuf2$emi_mar,noshuf2$emi_mar)),
     ylab=expression(paste("",F[ST],sep="")),bty="l",cex.axis=1.5,cex.lab=1.85,xlim=xr,xlab="",type="l",lwd=2,col="grey50",lty="dashed")
lines(con_m_fst~unique(con$t),lwd=2)
text(1030,par("usr")[[4]]*0.95,"C",cex=2)

mtext("time [generations]",side=1,line=.5,outer=T,cex=1.25,at=.575)
#title(xlab="time [generations]",line=.5,outer=T,cex.lab=1.85)
dev.off()