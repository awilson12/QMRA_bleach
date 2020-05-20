#Developed by Amanda M. Wilson and Mark H. Weir

#Bleach QMRA Code - relating bleach efficacies to estimated COVID-19 infection risks

#Exposure scenario: single hand-to-surface contact, followed by single hand-to-mouth contact (dose)

#set up directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

#install required packages if not already installed

if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}
if("ggplot2" %in% rownames(installed.packages())==FALSE){install.packages("ggplot2"); require(ggplot2)}else{require(ggplot2)}
if("ggpubr" %in% rownames(installed.packages())==FALSE){install.packages("ggpubr"); require(ggpubr)}else{require(ggpubr)}

#read in bootstrapped values for dose-response
exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')


COVIDqmra<-function(disinfection=c(TRUE),iter){
  
  if (disinfection==TRUE){
    #log10 reduction expected
    runif(iter,1,5) #explore a wide range for now to see how this affects infection risk
  }else{
    #this is a baseline scenario, so no disinfection, therefore 0 reduction
    reduce<-rep(0,iter)
  }
  
 concsurf<-10^runif(iter,-0.1,4) #place holder distribution for surface concentrations
 TE.SH<-runif(iter,0,1) #place holder transfer efficiencies (surface-->hand)
 SH<-runif(iter,0,1) #place holder fraction of total hand surface area used
 
 TE.HM<-runif(iter,0,1) #place holder transfer efficiencies (hand-->mouth)
 Ahand<-runif(iter,445,535) #hand surface areas (Beamer et al., 2015 office study and Exposure Factors Handbook)
 SH.mouth<-runif(iter,0.01,0.04) #place holder fraction of hand used for hand-to-mouth contact but approximately single 
                                 #fingertip to multiple fingertips
 
 #estimating virus transfer to hand for hand-to-surface contact
 conchand<-concsurf*TE*SH
 
 #estimating dose for hand-to-mouth contact
 dose<-conchand*TE.HM*Ahand*SH.mouth
 
 #initialize vector for storing infection risks
 
 for (i in 1:iter){
   pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
   infect[i]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[i], give=FALSE, strict=TRUE)
   
   if(infect[i]==0){
     infect[i]<-1*10^-15 #cannot have zero infection risk, so replace with small risk
     
    }
 }
 
 sim.frame<-data.frame(infect=infect,dose=dose,conchand=conchand,TE.HM=TE.HM,Ahand=Ahand,
                       SH.mouth=SH.mouth,SH=SH,TE.SH=TE.SH,concsurf=concsurf,)
  
}


