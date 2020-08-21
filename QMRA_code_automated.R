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

set.seed(34)

COVIDqmra<-function(disinfection=c(TRUE),iter,RNAinfective){
  
  if (disinfection==TRUE){
    #log10 reduction expected
    reduce<-runif(iter,1,5) #explore a wide range for now to see how this affects infection risk
  }else{
    #this is a baseline scenario, so no disinfection, therefore 0 reduction
    reduce<-rep(0,iter)
  }
  
  
  TE.SH<-runif(iter,0.01,0.406) 
  SH<-runif(iter,0.01,0.25) 
  
  TE.HM<-rtrunc(iter,"norm",a=0,b=1,mean=0.3390,sd=0.15)
  Ahand<-runif(iter,445,535) #hand surface areas (Beamer et al., 2015 office study and Exposure Factors Handbook)
  SH.mouth<-runif(iter,0.008,0.012) 
  
  A.fomite<-runif(iter,15,150)*929 #range from small to large room total amount of SA that desks account for, for example (in cm^2)
  #in ft^2 and converting to cm^2. Assuming 50% capacity of small classroom (10 seats) and large classroom (100 seats) with 3 ft^2 per person
  
  #concentration on surface in area we're contacting in gc/cm^2
  concsurfstart<-(10^runif(iter,-1,4))
  
  #adjustment after we account for fraction of our hand in contact, the fraction of virus that's infective,
  #the total area of fomites available for contact, and the log10 reduction if disinfection was used
  concsurf<-((concsurfstart*SH*Ahand*RNAinfective)/A.fomite)/(10^reduce)
  
  #concentration on our hand due to transfer efficiency
  conchand<-concsurf*TE.SH
 
  #estimating dose for hand-to-mouth contact
  dose<-conchand*TE.HM*Ahand*SH.mouth
 
  #initialize vector for storing infection risks
  infect<-rep(NA,iter)
  alpha<-rep(NA,iter)
  beta<-rep(NA,iter)
 
  #infection risk estimate with randomly sampled pair of alpha and beta
  for (i in 1:iter){
    #randomly sample 1 through the final position in the data frame
    pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
    #use this row number as the row for the alpha and beta pair to be used
    #arguments for teh hyperg_1F1 function are alpha, alpha + beta, and -dose
    infect[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[i], give=FALSE, strict=TRUE)
    
    #save alpha and beta used for later correlation coefficient 
    alpha[i]<-exactbp$alpha[pair]
    beta[i]<-exactbp$Beta[pair]
    
    #replace zero risks with small risk value
    if(infect[i]==0){
      infect[i]<-1*10^-15 #cannot have zero infection risk, so replace with small risk
    }
 }
 
 #save inputs and outputs for later spearman corr. coeff. 
 sim.frame<-data.frame(infect=infect,dose=dose,conchand=conchand,concsurfstart=concsurfstart,concsurf=concsurf,TE.HM=TE.HM,Ahand=Ahand,
                       SH.mouth=SH.mouth,SH=SH,TE.SH=TE.SH,RNAinfective=rep(RNAinfective,iter),
                       disinfect=rep(disinfection,iter),reduce=reduce,alpha=alpha,beta=beta,A.fomite=A.fomite)
 
 sim.frame<<-sim.frame
}

#-----------------automated run of sim function for scenarios-----------------------------------

SIM <- c("lowinfect_nodisinfect","highinfect_nodisinfect",
         "lowinfect_disinfect","highinfect_disinfect")   

NUM.SIM <- length(SIM)     

iter<-10000

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

View(cormat)

sim.frame.all$concenstatus<-rep(NA,length(sim.frame.all$infect))

#less than 1 gc/cm^2 in area we're contacting, low viral bioburden
#1 or more gc/cm^2 in area we're contacting, high viral bioburden
sim.frame.all$concenstatus[sim.frame.all$concsurfstart<1]<-"low"
sim.frame.all$concenstatus[sim.frame.all$concsurfstart>=1]<-"high"

sim.frame.all$RNAinfective[sim.frame.all$RNAinfective==0.01]<-"1% Infective"
sim.frame.all$RNAinfective[sim.frame.all$RNAinfective==0.1]<-"10% Infective"

sim.frame.all$reductionrange<-rep(NA,length(sim.frame.all$infect))
sim.frame.all$reductionrange[sim.frame.all$reduce==0]<-"0"
sim.frame.all$reductionrange[sim.frame.all$reduce>=1 & sim.frame.all$reduce<=2]<-"1-2"
sim.frame.all$reductionrange[sim.frame.all$reduce>2 & sim.frame.all$reduce<=3]<-"2-3"
sim.frame.all$reductionrange[sim.frame.all$reduce>3 & sim.frame.all$reduce<=4]<-"3-4"
sim.frame.all$reductionrange[sim.frame.all$reduce>4 & sim.frame.all$reduce<=5]<-"4-5"


#---- plots for brief note -------------------

windows()
ggplot(sim.frame.all)+geom_boxplot(aes(x=reductionrange,y=infect,group=interaction(concenstatus,reductionrange),fill=concenstatus))+
  facet_wrap(~RNAinfective)+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_x_discrete(name=expression("Log"[10]*phantom(x)*"Reduction"))+
  scale_fill_grey(name="Viral Bioburden",labels=c(expression("1 to 10,000 gc/cm"^2),expression("Less than 1 gc/cm"^2)),start=0.4,end=.8)+
  geom_hline(yintercept = 1e-4,linetype="dashed",size=1.5,colour="red")+
  geom_hline(yintercept = 1e-6,linetype="dashed",size=1.5,colour="orange")+
  theme_pubr()+theme_bw()+theme(axis.text=element_text(size=16),axis.title=element_text(size=16),
                                strip.text=element_text(size=16),
                                legend.text=element_text(size=16),
                                axis.text.x = element_text(angle = 45),
                                legend.title=element_text(size=16))

windows()
ggplot(sim.frame.all)+geom_boxplot(aes(x=concenstatus,y=infect,group=interaction(concenstatus,disinfect),fill=disinfect))+
  facet_wrap(~RNAinfective)+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_x_discrete(name="Bioburden",labels=c(expression("1 to 10,000 gc/cm"^2),expression("Less than 1 gc/cm"^2)))+
  scale_fill_grey(name="",labels=c(expression("No log"[10]*phantom(x)*"reduction"),expression("1-5 log"[10]*phantom(x)*"reduction")),start=0.4,end=.8)+
  geom_hline(yintercept = 1e-4,linetype="dashed",size=1.5,colour="red")+
  geom_hline(yintercept = 1e-6,linetype="dashed",size=1.5,colour="orange")+
  theme_pubr()+theme_bw()+theme(axis.text=element_text(size=16),axis.title=element_text(size=16),
                                strip.text=element_text(size=16),
                                legend.text=element_text(size=16),
                                legend.title=element_text(size=16))