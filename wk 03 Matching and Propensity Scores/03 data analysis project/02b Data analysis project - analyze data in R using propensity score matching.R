# For this assignment we will use data from Lalonde (1986),
# that aimed to evaluate the impact of  National
# Supported Work (NSW) Demonstration, which is a labor training program, on
# post-intervention income levels. Interest is in estimating the causal effect of
# this training program on income.

# load packages

#install.packages(“tableone”)
#install.packages(“Matching”)
#install.packages("MatchIt")

library(tableone)
library(Matching)
library(MatchIt)


#import data
data(lalonde, package = 'MatchIt')



#create variables
age<-lalonde$age
educ<-lalonde$educ
black<-as.numeric(lalonde$race=='black')
hispan<-as.numeric(lalonde$race=='hispan')
married<-lalonde$married
nodegr<-lalonde$nodegr
re74<-lalonde$re74
re75<-lalonde$re75
re78<-lalonde$re78
treat<-lalonde$treat

mydata=cbind(age,educ,black,hispan,married,nodegr,re74,re75,treat,re78)
mydata<-data.frame(mydata)

#covariates we will use
xvars<-c("age","educ","black","hispan","married","nodegr","re74","re75")

#look at a table 1 - standard differences for all of the confounding variables (pre-matching)
table1<- CreateTableOne(vars=xvars,strata="treat", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)


#fit a propensity score model. logistic regression
psmodel<-glm(treat ~ age + educ + black + hispan + married + nodegr + re74 + re75,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


# Now carry out propensity score matching using the Match function. 
# Match on the propensity score itself, not logit of the propensity score.  Obtain the standardized differences for the matched data


# result 1 - without caliper
set.seed(931139)
logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treat,M=1,X=pscore,replace=FALSE) #without caliper
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("age","educ","black","hispan","married","nodegr","re74","re75","re78")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treat", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

# result 2 - with caliper
set.seed(931139)
logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treat,M=1,X=pscore,replace=FALSE, caliper=0.1) #with caliper
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("age","educ","black","hispan","married","nodegr","re74","re75","re78")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treat", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$re78[matched$treat==1]
y_con<-matched$re78[matched$treat==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)
