#Developed by Amanda M. Wilson and Mark H. Weir

#Bleach QMRA Code - relating bleach efficacies to estimated COVID-19 infection risks

#Exposure scenario: single hand-to-surface contact, followed by single hand-to-mouth contact (dose)

#set up directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

#install required packages if not already installed
if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}
if("ggplot2" %in% rownames(installed.packages())==FALSE){install.packages("ggplot2"); require(ggplot2)}else{require(ggplot2)}
if("ggpubr" %in% rownames(installed.packages())==FALSE){install.packages("ggpubr"); require(ggpubr)}else{require(ggpubr)}
if("reshape2" %in% rownames(installed.packages())==FALSE){install.packages("reshape2"); require(reshape2)}else{require(reshape2)}


#read in bootstrapped values for dose-response
exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')


COVIDqmra<-function(disinfection=c(TRUE),iter,RNAinfective){
  
  if (disinfection==TRUE){
    #log10 reduction expected
    reduce<-runif(iter,1,5) #explore a wide range for now to see how this affects infection risk
  }else{
    #this is a baseline scenario, so no disinfection, therefore 0 reduction
    reduce<-rep(0,iter)
  }
  
 concsurf<-10^runif(iter,-0.1,4)*RNAinfective #place holder distribution for surface concentrations
 TE.SH<-runif(iter,0.01,0.406) #place holder transfer efficiencies (surface-->hand)
 SH<-runif(iter,0.01,0.25) #place holder fraction of total hand surface area used
 
 TE.HM<-rtrunc(iter,"norm",a=0,b=1,mean=0.3390,sd=0.15)
 Ahand<-runif(iter,445,535) #hand surface areas (Beamer et al., 2015 office study and Exposure Factors Handbook)
 SH.mouth<-runif(iter,0.008,0.012) #place holder fraction of hand used for hand-to-mouth contact but approximately single 
                                 #fingertip to multiple fingertips
 
 #estimating virus transfer to hand for hand-to-surface contact
 conchand<-(concsurf*TE.SH*SH)/10^reduce
 
 #estimating dose for hand-to-mouth contact
 dose<-conchand*TE.HM*Ahand*SH.mouth
 
 #initialize vector for storing infection risks
 infect<-rep(NA,iter)
 alpha<-rep(NA,iter)
 beta<-rep(NA,iter)
 
 for (i in 1:iter){
   pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
   infect[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[i], give=FALSE, strict=TRUE)
   alpha[i]<-exactbp$alpha[pair]
   beta[i]<-exactbp$Beta[pair]
     
   if(infect[i]==0){
     infect[i]<-1*10^-15 #cannot have zero infection risk, so replace with small risk
    }
 }
 
 sim.frame<-data.frame(infect=infect,dose=dose,conchand=conchand,TE.HM=TE.HM,Ahand=Ahand,
                       SH.mouth=SH.mouth,SH=SH,TE.SH=TE.SH,concsurf=concsurf,RNAinfective=rep(RNAinfective,iter),
                       disinfect=rep(disinfection,iter),reduce=reduce,alpha=alpha,beta=beta)
 
 sim.frame<<-sim.frame
}

#-----------------automated run of sim function for scenarios-----------------------------------

SIM <- c("lowinfect_nodisinfect","highinfect_nodisinfect",
         "lowinfect_disinfect","highinfect_disinfect")   

NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

iter<-1000

for(j in 1:NUM.SIM)
  
{
  sim.num <- j; sim.name <- SIM[j]
  
  if(sim.name=="lowinfect_nodisinfect"){
    if(dir.exists("lowinfect_nodisinfect")==FALSE)
    {dir.create("lowinfect_nodisinfect"); setwd("lowinfect_nodisinfect")}
    else{setwd("lowinfect_nodisinfect")}
    
    RNAinfect<-0.01
    disinfection<-c(FALSE)
    
    }
  
  if(sim.name=="highinfect_nodisinfect"){
    if(dir.exists("highinfect_nodisinfect")==FALSE)
    {dir.create("highinfect_nodisinfect"); setwd("highinfect_nodisinfect")}
    else{setwd("highinfect_nodisinfect")}
    
    RNAinfect<-0.1
    disinfection<-c(FALSE)
    
    }
  
  if(sim.name=="lowinfect_disinfect"){
    if(dir.exists("lowinfect_disinfect")==FALSE)
    {dir.create("lowinfect_disinfect"); setwd("lowinfect_disinfect")}
    else{setwd("lowinfect_disinfect")}
    
    RNAinfect<-0.01
    disinfection<-c(TRUE)
    
    }
  if(sim.name=="highinfect_disinfect"){
    if(dir.exists("highinfect_disinfect")==FALSE)
    {dir.create("highinfect_disinfect"); setwd("highinfect_disinfect")}
    else{setwd("highinfect_disinfect")}
    
    RNAinfect<-0.1
    disinfection<-c(TRUE)
    
    }
  
  COVIDqmra(disinfection=disinfection,RNAinfective = RNAinfect,iter=iter)
  
  write.csv(sim.frame,file=sprintf("%s.sim.frame.csv",sim.name))
  
  #reset directory to parent folder so we can go to correct subfolder within parent folder for next sim run
  setwd(this.dir)
  
  if (j==1){
   sim.frame.all<-sim.frame
  }else{
    sim.frame.all<-rbind(sim.frame.all,sim.frame)
  }
  
}


framecor = sim.frame.all

cor(sim.frame.all$infect,sim.frame.all$Ahand,method=c("spearman"))

cormat<-cor(framecor,method=c("spearman"))
melted_cormat<-melt(cormat)
ggplot(data=melted_cormat,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
  geom_text(aes(label = signif(value, 2))) +
  scale_fill_gradient(low = "white", high = "blue") 


sim.frame.all$concenstatus<-rep(NA,length(sim.frame.all$infect))
sim.frame.all$concenstatus[sim.frame.all$concsurf<1]<-"low"
sim.frame.all$concenstatus[sim.frame.all$concsurf>=1]<-"high"

sim.frame.all$RNAinfective[sim.frame.all$RNAinfective==0.01]<-"1% Infective"
sim.frame.all$RNAinfective[sim.frame.all$RNAinfective==0.1]<-"10% Infective"


sim.frame.all$reductionrange<-rep(NA,length(sim.frame.all$infect))
sim.frame.all$reductionrange[sim.frame.all$reduce==0]<-"0"
sim.frame.all$reductionrange[sim.frame.all$reduce>0 & sim.frame.all$reduce<=1]<-"0-1"
sim.frame.all$reductionrange[sim.frame.all$reduce>1 & sim.frame.all$reduce<=2]<-"1-2"
sim.frame.all$reductionrange[sim.frame.all$reduce>2 & sim.frame.all$reduce<=3]<-"2-3"
sim.frame.all$reductionrange[sim.frame.all$reduce>3 & sim.frame.all$reduce<=4]<-"3-4"
sim.frame.all$reductionrange[sim.frame.all$reduce>4 & sim.frame.all$reduce<=5]<-"4-5"


#---- plots for brief note -------------------

windows()
ggplot(sim.frame.all)+geom_boxplot(aes(x=reductionrange,y=infect,group=interaction(concenstatus,reductionrange),fill=concenstatus))+
  facet_wrap(~RNAinfective)+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_x_discrete(name=expression("Log"[10]*phantom(x)*"Reduction"))+
  scale_fill_grey(name="Contamination",labels=c(expression(phantom(x)>="1 RNA/cm"^2),expression("< 1 RNA/cm"^2)),start=0.4,end=.8)+
  geom_hline(yintercept = 1e-4,linetype="dashed",size=1,colour="red")+
  geom_hline(yintercept = 1e-6,linetype="dashed",size=1,colour="orange")+
  theme_pubr()+theme_bw()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
                                strip.text=element_text(size=12),
                                legend.text=element_text(size=12),
                                axis.text.x = element_text(angle = 45))

windows()
ggplot(sim.frame.all)+geom_boxplot(aes(x=concenstatus,y=infect,group=interaction(concenstatus,disinfect),fill=disinfect))+
  facet_wrap(~RNAinfective)+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_x_discrete(name="Contamination",labels=c(expression(phantom(x)>="1 RNA/cm"^2),expression("< 1 RNA/cm"^2)))+
  scale_fill_grey(name="",labels=c("No disinfection","Disinfection"),start=0.4,end=.8)+
  geom_hline(yintercept = 1e-4,linetype="dashed",size=1,colour="red")+
  geom_hline(yintercept = 1e-6,linetype="dashed",size=1,colour="orange")+
  theme_pubr()+theme_bw()+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
                                strip.text=element_text(size=12),
                                legend.text=element_text(size=12))

#---- exploratory plots-------------------

ggplot(sim.frame.all)+stat_ecdf(aes(x=infect,group=interaction(disinfect,concenstatus),colour=disinfect,linetype=concenstatus),size=1)+
  scale_x_continuous(trans="log10",name="Infection Risk")+
  scale_colour_discrete(name="",labels=c("Surfaces Not Disinfected","Surfaces Disinfected"))+
  scale_linetype_discrete(name="",labels=c(expression(">= 1 viral particle/cm"^2),expression("<1 viral particle/cm"^2)))+
  scale_y_continuous(name="Fraction of Data")+
  theme_pubr()+theme(legend.position = "right")
#geom_vline(xintercept=1e-4)+
#geom_vline(xintercept=1e-6)

A<-ggplot(sim.frame.all[sim.frame.all$reduce>0,])+geom_point(aes(x=concsurf,y=infect,colour=reduce))+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_colour_continuous(name=expression("log"[10]*phantom(x)*"reduction"))+
  scale_x_continuous(trans="log10",name=expression("Surface Concentration (viral particles/cm"^2*")"))+
  geom_hline(yintercept=1e-4,linetype="dashed",colour="red",size=1.5)+
  theme_pubr()+ggtitle("1/10,000 Risk Target")

B<-ggplot(sim.frame.all[sim.frame.all$reduce>0,])+geom_point(aes(x=concsurf,y=infect,colour=reduce))+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_colour_continuous(name=expression("log"[10]*phantom(x)*"reduction"))+
  scale_x_continuous(trans="log10",name=expression("Surface Concentration (viral particles/cm"^2*")"))+
  geom_hline(yintercept=1e-6,linetype="dashed",colour="red",size=1.5)+
  theme_pubr()+ggtitle("1/1,000,000 Risk Target")

ggarrange(A,B,common.legend = TRUE)
