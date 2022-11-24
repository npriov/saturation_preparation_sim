library(MASS)
library(cowplot)
library(RColorBrewer)
library(pracma)
library(ggplot2)

throt<- function(phi,theta){
  Rz = zrot(-theta);
  Rx = xrot(phi);
  Rth = ginv(Rz) %*% Rx %*% Rz;
  return(Rth);}
xrot<- function(phi){
  Rx = t(matrix(c(1,0,0,0, cos(phi), -sin(phi),0, sin(phi), cos(phi)),nrow=3,ncol=3));
  return(Rx);}
zrot<- function(phi){
Rz = t(matrix(c(cos(phi), -sin(phi), 0,sin(phi), cos(phi), 0,0, 0, 1),nrow=3,ncol=3));
return(Rz);}
freeprecess <- function(T,T1,T2,df){
  
  #	Function simulates free precession and decay
  #	over a time interval T, given relaxation times T1 and T2
  #	and off-resonance df.  Times in ms, off-resonance in Hz.
  
  phi = 2*pi*df*T/1000;	# Resonant precession, radians.
  E1 = exp(-T/T1);	
  E2 = exp(-T/T2);
  
  Afp=t(matrix(c(E2, 0, 0,0, E2, 0, 0, 0, E1),nrow=3,ncol=3))%*%zrot(phi);
  
  Bfp = c(0, 0, 1-E1);

  return(list(Afp,Bfp));
}
sat_effect <- function(rois_cases){
  # Function that simulates the effect of the saturation pulse
  
  #define saturation train
  FA=155; dT=0.005; #in ms
  dur=0.5; sub=rep(1,length(seq(0,dur,dT))-1);
  b1=c(-sub,sub,sub,-sub,-sub,sub,sub,-sub,-sub,sub,sub,-sub); 
  rf=b1*FA*pi*dT/180;
  freq=0;
  
  #rotations
  M=rois_cases[[2]];
  M = throt(abs(rf[1]),angle(rf[1])) %*% M;	# RF Rotation.
  for (rf_inst in 2:length(rf)){
    AGH= freeprecess(dT,rois_cases[[3]],rois_cases[[4]],0);
    M = AGH[[1]]%*%M+AGH[[2]];				# Propagate to next pulse.
    M = throt(abs(rf[rf_inst]),angle(rf[rf_inst])) %*% M;	# RF Rotation.
  }
  rois_cases[[2]]=M;
  return(rois_cases)
}



n_font=10

pulsetrains<-1:300
sat_to_sat=0.381;

# ROI, f, k, Rm, Rw, Sm0,  T2 of free water in ms 
k=2; Rm=2; Rw=0.4; Sm0=1-0.075; 
MT_init_values<-list(list("WM",0.30, k, Rm, Rw, Sm0,900),
                     list("GM",0.15, k, Rm, Rw, Sm0,900),
                     list("AR",0.001, k, Rm, Rw, Sm0,67.5),
                     list("VE",0.001, k, Rm, Rw, Sm0,7),
                     list("CSF",0.001, k, Rm, Rw, Sm0,900))



simulationRm_df_ts<-list()
simulationRm_df_ts_betweentrains<-list()

t=seq(0.001,sat_to_sat,length.out=15)

for (roi in 1:length(MT_init_values)){
  M0_start=c(0,0,1); #initialize M0
  out_lit<-sat_effect(list("FW",M0_start,(1/MT_init_values[[roi]][[5]])*1000,MT_init_values[[roi]][[7]])) #run first saturation pulse on FW
  Sw0=1-out_lit[[2]][3]; #transform to FW fractional saturation from Mz/M0 ratio
  Sm0=MT_init_values[[roi]][[6]]; #grab MW fractional saturation as well
  f1=MT_init_values[[roi]][[2]]; #initialize bound pool size
  Rw=MT_init_values[[roi]][[4]]; #initiWlize R1 of FW
  Rm=MT_init_values[[roi]][[5]]; #initialize R1 of MM
  
  km=k/f1; #calculate bound to water exchange rate
  kw=k/(1-f1); #calculate water to bound exchange rate
  l1=(Rw+Rm+kw+km+sqrt((Rm-Rw+km-kw)^2+4*kw*km))/2; #cross-over relaxation rate 1
  l2=(Rw+Rm+kw+km-sqrt((Rm-Rw+km-kw)^2+4*kw*km))/2; #cross-over relaxation rate 2

  Mzbetweentrains<-rep(0,length(pulsetrains))
  for (n in pulsetrains){  
    a1=((Rw+kw-l2)*Sw0-kw*Sm0)/(l1-l2); #update the amplitudes based on current saturation level
    a2=-((Rw+kw-l1)*Sw0-kw*Sm0)/(l1-l2);
    y1 =a1*exp(-l1*t)+a2*exp(-l2*t); #free MT exchange after saturation
    Mzbetweentrains[n]<-1-y1[5]; #save the initial Mz for each pulse train
    
    M0_start[3]=tail(1-y1,1); #reinitialize Mz/M0 for new saturation train
    out_lit<-sat_effect(list("FW",M0_start,(1/MT_init_values[[roi]][[5]])*1000,MT_init_values[[roi]][[7]])); #apply new saturation pulse 
    M0_start<-as.vector(out_lit[[2]]);
    Sw0=1-M0_start[3];
    
  }
  simulationRm_df_ts_betweentrains[[roi]]<-data.frame(Mzbetweentrains,pulsetrains,MT_init_values[[roi]][[1]])
  simulationRm_df_ts[[roi]]<-data.frame(1-y1,t,Rw,f1,MT_init_values[[roi]][[1]])
}



simulationRm_df_ts_betweentrains<-do.call(rbind, simulationRm_df_ts_betweentrains)
colnames(simulationRm_df_ts_betweentrains)<-c("Mz","n","ROI")
simulationRm_df_ts_betweentrains$ROI<-as.factor(simulationRm_df_ts_betweentrains$ROI)

simulationRm_df_ts_betweentrains <- subset(simulationRm_df_ts_betweentrains, n<=20)






simRw_trains<-ggplotGrob(ggplot(data=simulationRm_df_ts_betweentrains,aes(x=n,y=Mz,color=as.factor(ROI),group=as.factor(ROI)))+
                    geom_line(aes(color=ROI),size=1)+
                    theme_classic()+ylab(expression(M["z"]*" free water (ratio)"))+xlab("Saturation trains")+
                    scale_color_manual(labels = c("WM", "GM","AR","VE","CSF"), values = c("#009E73", "#E69F00","#AE0000","#007ab899","#000000")) +
                    ylim(0.4, 1)+
                    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),text = element_text(size=8))+
                    theme(axis.text=element_text(size=n_font),axis.title=element_text(size=n_font))+  #guides(color=FALSE)+
                    theme(legend.title=element_text(size=n_font),legend.direction="vertical",legend.text=element_text(size=n_font))+theme(legend.box = "vertical")+
                    theme(legend.title = element_blank())+guides(color=FALSE))
ggdraw(simRw_trains)



simulationRm_df<-simulationRm_df_ts
simulationRm_df<-do.call(rbind, simulationRm_df)  
colnames(simulationRm_df)<-c("y1","t","Tbound","f","ROI")
simulationRm_df$Tbound<-as.factor(simulationRm_df$Tbound)
simulationRm_df$ROI<-as.factor(simulationRm_df$ROI)

simRw<-ggplotGrob(ggplot(data=simulationRm_df,aes(x=t*1000,y=y1,color=as.factor(ROI),group=as.factor(ROI)))+
                    geom_line(aes(color=ROI),size=1)+
                    theme_classic()+ylab("Mz free water (ratio)")+xlab("Time between saturations (ms)")+
                    scale_color_manual(labels = c("WM", "GM","AR","VE","CSF"), values = c("#009E73", "#E69F00","#AE0000","#007ab899","#000000")) +
                    ylim(0.4, 1)+
                    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),text = element_text(size=8))+
                    theme(axis.text=element_text(size=n_font),axis.title=element_text(size=n_font))+  #guides(color=FALSE)+
                    theme(legend.title=element_text(size=n_font),legend.direction="vertical",legend.text=element_text(size=n_font))+theme(legend.box = "vertical")+
                    theme(legend.title = element_blank())+guides(color=FALSE))
ggdraw(simRw)

