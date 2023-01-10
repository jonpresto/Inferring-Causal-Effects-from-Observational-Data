###################
#RHC Example
#https://biostat.app.vumc.org/wiki/pub/Main/DataSets/rhc.html 

#install packages
#install.packages("tableone")
#devtools::install_github(repo = "kaz-yos/tableone", ref = "develop")
#install.packages("Matching")
#install.packages("MatchIt")

#load packages
library(tableone)
library(MatchIt)
library(Matching)

#read in data
#load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
load(url("https://biostat.app.vumc.org/wiki/pub/Main/DataSets/rhc.sav"))
#view data
View(rhc)

#treatment variables is swang1
#x variables that we will use
#cat1: primary disease category
#age
#sex
#meanbp1: mean blood pressure

#create a data set with just these variables, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
died<-as.numeric(rhc$death=='Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

#new dataset
mydata<-cbind(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,
              age,female,meanbp1,treatment,died)
mydata<-data.frame(mydata)

#covariates we will use (shorter list than you would use in practice)
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

########################################################
#option 1: direct way of creating propensity score
########################################################
#fit a propensity score model. logistic regression
psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values

#combine pscore with mydata
mydata2 = cbind(mydata,pscore)

library(ggplot2)
ggplot(mydata2, aes(x=pscore)) +
  geom_histogram(fill="white", colour="black") +
  facet_grid(treatment ~ .)

########################################################
#option 2: use MatchIt package propensity score
########################################################
#use matchit for propensity score, nearest neighbor matching
m.out <- matchit(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+sepsis+age+female+meanbp1,
                 data=mydata, method="nearest")

summary(m.out)

#propensity score plots
plot(m.out,type="jitter")
plot(m.out,type="hist")


#look at a table 1
table1<- CreateTableOne(vars=xvars,strata="treatment", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)

############################################
#do greedy matching on logit(PS) using Match without a caliper
############################################
logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)


############################################
#do greedy matching on logit(PS) using Match with a caliper
############################################
logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#get standardized differences
matchedtab2<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab2, smd = TRUE)

