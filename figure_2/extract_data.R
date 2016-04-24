exps<-c("K_150","K_50","lambda_15","lambda_3","mu_015","mu_03","mu_045","sigma_075","sigma_025")
emi_mar_shuf_m<-numeric()
emi_mar_shuf_q025<-numeric()
emi_mar_shuf_q075<-numeric()
emi_mar_noshuf_m<-numeric()
emi_mar_noshuf_q025<-numeric()
emi_mar_noshuf_q075<-numeric()
velo_shuf_m<-numeric()
velo_shuf_q025<-numeric()
velo_shuf_q075<-numeric()
velo_noshuf_m<-numeric()
velo_noshuf_sd<-numeric()
velo_noshuf_q025<-numeric()
velo_noshuf_q075<-numeric()

velo_diff_m<-numeric()
velo_diff_q025<-numeric()
velo_diff_q075<-numeric()
emi_diff_m<-numeric()
emi_diff_q025<-numeric()
emi_diff_q075<-numeric()

for (i in 1:length(exps)) {
  con <- read.table(paste(exps[i],"/control/data/output_noshuffle.txt",sep=""),header=T)
  shuf <-read.table(paste(exps[i],"/shuffle/data/output_shuffle.txt"  ,sep=""),header=T)
  
  emis_mar_shuf<-numeric()
  emis_mar_noshuf<-numeric()
  velos_shuf<-numeric()
  velos_noshuf<-numeric()
  velos_diff<-numeric()
  emis_diff<-numeric()
  
  for (r in unique(con$repli)) {
    
    maxrange <- max(con$real_b[con$repli==r])
    fill_time <- con$t[which((con$repli==r)&(con$real_b==maxrange))][[1]]
    
    velos_noshuf[r] <- con$real_b[which((con$t==fill_time)&(con$repli==r))]/(fill_time-1000)
    velos_shuf[r] <- shuf$real_b[which((shuf$t==fill_time)&(shuf$repli==r))]/(fill_time-1000)
    
    velos_diff[r] <- (velos_noshuf[r]-velos_shuf[r])/(velos_noshuf[r])
    
    emis_mar_noshuf[r] <- con$emi_mar[which((con$t==fill_time)&(con$repli==r))]
    emis_mar_shuf[r] <- shuf$emi_mar[which((shuf$t==fill_time)&(shuf$repli==r))]
    
    emis_diff[r] <- (emis_mar_noshuf[r]-emis_mar_shuf[r])/(emis_mar_noshuf[r])
    
  }
  
  velo_shuf_m[i] <- median(velos_shuf,na.rm=T)
  velo_shuf_q025[i] <- quantile(velos_shuf,probs=0.25,na.rm=T)
  velo_shuf_q075[i] <- quantile(velos_shuf,probs=0.75,na.rm=T)
  velo_noshuf_m[i] <- median(velos_noshuf,na.rm=T)
  velo_noshuf_q025[i] <- quantile(velos_noshuf,probs=0.25,na.rm=T)
  velo_noshuf_q075[i] <- quantile(velos_noshuf,probs=0.75,na.rm=T)
  
  emi_mar_shuf_m[i] <- median(emis_mar_shuf,na.rm=T)
  emi_mar_shuf_q025[i] <- quantile(emis_mar_shuf,probs=0.25,na.rm=T)
  emi_mar_shuf_q075[i] <- quantile(emis_mar_shuf,probs=0.75,na.rm=T)
  emi_mar_noshuf_m[i] <- median(emis_mar_noshuf,na.rm=T)
  emi_mar_noshuf_q025[i] <- quantile(emis_mar_noshuf,probs=0.25,na.rm=T)
  emi_mar_noshuf_q075[i] <- quantile(emis_mar_noshuf,probs=0.75,na.rm=T)
  
  velo_diff_m[i] <- median(velos_diff,na.rm=T)
  velo_diff_q025[i] <- quantile(velos_diff,probs=0.25,na.rm=T)
  velo_diff_q075[i] <- quantile(velos_diff,probs=0.75,na.rm=T)
  emi_diff_m[i] <- median(emis_diff,na.rm=T)
  emi_diff_q025[i] <- quantile(emis_diff,probs=0.25,na.rm=T)
  emi_diff_q075[i] <- quantile(emis_diff,probs=0.75,na.rm=T)
}

pars<-c("K","K","lambda","lambda","mu","mu","mu","sigma","sigma")
vals<-c(150,50,1.5,3,0.15,0.3,0.45,0.75,0.25)

data<-data.frame(exps,pars,vals,velo_shuf_m,velo_shuf_q025,velo_shuf_q075,velo_noshuf_m,velo_noshuf_q025,velo_noshuf_q075,emi_mar_shuf_m,emi_mar_shuf_q025,emi_mar_shuf_q075,emi_mar_noshuf_m,emi_mar_noshuf_q025,emi_mar_noshuf_q075,velo_diff_m,velo_diff_q025,velo_diff_q075,emi_diff_m,emi_diff_q025,emi_diff_q075)

write.table(data,file="results.dat",quote=F)