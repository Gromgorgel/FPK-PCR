## Disclaimer: This function has only been tested on a limited amount of data
##             there may be a few bugs left. Please contact the author with 

##             any remarks and/or suggestions.
##
## email: antoon.lievens@ec.europa.eu
##
## Most problems are caused by a lack of data resulting in a suboptimal fit
## try running more cycles and increasing primer concentration.
## If that doesn't help, feel free to contact me.

################################################################################
#The V13 series is based on V6-11 and searches for a better Emax and alpha estimates
#Version 13 is the final outcome of the 11&12 series of experiments it has:
#     *  the baseline optimisation system
#     *  a limited freedom fit instead of regular quadratic
#     *  NO weighing or badness system for the  baseline optimisation
#     *  Kalman covariance via bootstrap 
#     *  resec implementation to weed out non-amplification fits
#     *  stada no longer cycle(5%)+1 but now cycle(7.5%)
#     *  debugg option added 
################################################################################

##LEGENDA
##plots = TRUE to make some nice output plots which may help locate problems
##baseline =  either "optimized", "slanted", or "flat" depending on the model to be used
##output =  either "estimates" (for E, ai0 & Cq), "parameters" (for model parameters as a list), or "all" (both as a concatenated vector)
##silent = TRUE to prevent the algorithm from printing informative messages
##kalman = TRUE to use kalman filter and include groundphase data into the analysis
##debugg = TRUE to report intermediate values, warning messages, and diagnosic plots
##resec.tolerance = 5PLM residual error per cycle cut-off, serves as a tool in low-level presence analysis (reactions with high values are not likely to contain bona fide amplification, NA is reported for reactions that exceed the threshold)
##BOOTSTRAPPING
##bootstrap = TRUE to enable bootstrap
##n = number of bootstraps
##p = p-value of the confidence bounds calculated

analyse<-function(fluo,base.line="optimized",output="estimates",bootstrap=FALSE,plots=FALSE,silent=FALSE,kalman=TRUE,debugg=FALSE,n=1000,p=0.05,resec.tolerance=0.125){ 
    fluo<-as.numeric(fluo)                  #coerce numeric format
    endd<-length(fluo)                      #number of cycles ran
    cycles<-c(1:endd)                       #generate cycle numbers      
    plotto<-1                               #indicates which type of plots are generated (depends on available data, default=1)
    fail<-0                                 #keeps track of failures in order to return appropriate error message (default=0)
    vfail<-0
    mark<-FALSE                             #keeps track of over background subtraction
    require(minpack.lm)                     #load Levenberg-Marquardt nonlinear least-squares algorithm 
    require(GenSA)
    require(MASS)

#print a long line (this enables to distinguish between messages from different analysis
#in a multi-reaction analysis (eg: apply(fluo,2,analyse))
if(!isTRUE(silent)){message("------------")}

################################################################################     
                ########################################
                ##         Start of Procedure         ##
                ##      Initial Control Section       ##
                ########################################
################################################################################     

########################################
##            Syntax-control          ##              
########################################

#baseline Syntax: check if user input is correct (model is either "slanted" or "flat")
if(any(base.line==c("slanted","flat","optimized"))){
#Output Syntax: check if user input is correct (output is either "estimates" , "parameters", or "all")
if(any(output==c("estimates","parameters","all"))){
#Logical Syntax: check if user input is correct (i.e. are the logicals indeed logicals)
if(all(is.logical(c(plots,silent,debugg)))){

#for debugging: supress standard plots, allow info messages
if(isTRUE(debugg)&isTRUE(plots)){plots<-FALSE;silent<-FALSE}

########################################
##            Data-Control            ##
########################################

#check if there is an obvious absence of amplification (i.e. net increase in 
#fluorescence is less than twofold, or mean fluorescence is negative)

#note that we use a rough baseline estimate & remove all NAs (we cannot use endd since there might be trailing NAs!)
if(!((2*(mean(head(na.omit(fluo),n=5))-mean(head(na.omit(fluo),n=3))))>(mean(tail(na.omit(fluo),n=5))-mean(head(na.omit(fluo),n=3))))&!mean(fluo,na.rm=T)<0){ 


################################################################################     
                ########################################
                ##        End Control Section         ##
                ##    Start of Artifact correction    ##
                ########################################
################################################################################     
        

#notify user that data artifacts are being removed
if(min(fluo)<0|which.max(fluo)!=endd|any(is.na(fluo))){
    if(!isTRUE(silent)){message("Data artifacts found. Carying out pre-analysis data treatment:")}

########################################
##            Missing values          ##
########################################

if(any(is.na(fluo))){
  #Warn User
  if(!isTRUE(silent)){message("Data contains NA values...")}
  #check if last value is NA
  if(is.na(fluo[endd])){
    #check if all NA are consecutive
    if(all(diff(c(1:endd)[is.na(fluo)])==1)){
      #if so we can simply remove all NA from the dataset
      #trailing NAs do not affect analysis process
      fluo<-as.numeric(na.omit(fluo))
      #inform user
      if(!isTRUE(silent)){message("...Trailing NA values removed")}
    }else{
      #the most complicated case: trailing NAs + missing values
      #we only remove the trailing NAs
      fluo<-fluo[1:(endd-(length(fluo[is.na(fluo)][diff(c(1:endd)[is.na(fluo)])==1])+1))]
      #inform user
      if(!isTRUE(silent)){message("...Trailing NA values removed");message("...additional missing values detected: some cycles have NA fluorescence");message("...no action taken")}
      }
    #Since NAs have been removed we need to update some values
    endd<-length(fluo)                      #number of cycles ran
    cycles<-c(1:endd)                      #generate cycle numbers
    #if the last value is not an NA we have missing values
    }else{
    if(!isTRUE(silent)){message("...Missing values detected: some cycles have NA fluorescence");message("...no action taken")}
    }}
    
########################################
##    Over-background-subtraction     ##
########################################

#background subtraction (as carried out  by some hardware platforms)
#can result in negative fluorescence values. This poses a severe problem
#when evaluating the baseline. Moreover, negative fluorescence values
#are complete nonsense. There is no black hole in the reaction.
#Therefore addition of a constant (dummy fluorescence) to every fluorescence 
##value in order to undo this 'over-background-subtraction' is justified.

if(min(fluo,na.rm=T)<0){
  if(!isTRUE(silent)){message("Correcting for over-background-subtraction...")}
  dummy<-1+ceiling(abs(min(fluo,na.rm=T)))   
  fluo<-fluo+dummy
  mark<-TRUE  #create marker to revert changes before output
  }else{mark<-FALSE} 

########################################
##           Plateau Decline          ##
########################################

##fluorescence decline in the plateau is problematic for the efficiency model
##since this means that efficiency will rise again. Therefor, if decline is 
##found, the plateau is flattened

if(which.max(fluo)!=endd){
  if(!isTRUE(silent)){message("Correcting for plateau decline in raw data...")}
  frank<-which.max(fluo)
  plat<-max(fluo)+0.0001*c(1:(endd-frank))
  fluo<-c(fluo[1:frank],plat)
  }      

#mark end of pre-analysis steps
if(!isTRUE(silent)&any(exists("frank"),exists("dummy"))){message("----")}}

########################################
##        Rescaling of the data       ##
########################################

##this step has no effect whatsoever on the analysis process
##the only reason it is carried out is to prevent "system is computationally singular" error 
##when solving the X'X matrix during the Kalman filter. This may happen if the fluorescence reaches very high values (~10000 FU).
##because of the square in the quadratic equation, these values dominate the matrix and cause the error. Hence rescaling

if(min(fluo)<1){rudy<-10^-(floor(log10(min(fluo))))}else{rudy<-1} #prevent fluorescence smaller than 1 (this will screw with our second rescaling)
fluo<-rudy*fluo                                                   #first rescaling
rudy<-append(rudy,c(floor(min(fluo)),(ceiling(max(fluo))-floor(min(fluo)))/10)) #we store all values to be able to-denormalize everyting
fluo<-sapply(fluo,function(x){return((x-rudy[2])/rudy[3])})        #all values now run between 0 and 10, but don't start exactly at 0 and 10

################################################################################     
                ########################################
                ##     End of Artifact correction     ##
                ##    Start of Pre-Analysis Section   ##
                ########################################
################################################################################     

########################################
##     fit of 5PLM to raw data        ##
########################################
 
##get initial values (we use initial values from 4PLM and use g = 1 and s = 0.1)
  sp<-try(getInitial(fluo~SSfpl(cycles,a,d,xmid,b),data=data.frame(cycles,fluo)),silent=TRUE)
  if(inherits(sp, "try-error")){
    #failure of the initial value routine may happen in very "late" reactions (Cq>40)
    #However, this does not mean the model cannot be fitted
    #we fall back on manual calculation of the initial parameter values & inform user
    if(!isTRUE(silent)){message("getInitial(~SSfpl) failed");message("Switching to manual calculation...")}  
  	a<-max(fluo)               				  #calcute start parameter a
  	Y0<-mean(fluo[3:8])					        #calcute start parameter Y0
	  X0<-length(fluo[fluo<(Y0+(a-Y0)/2)])#calcute start parameter X0
    ##special treatment for b
    bgrid <- expand.grid(Y0= Y0, a= a, X0= X0, b= seq(0,20,by=0.5))		           #define the grid to search best b in
    FPLM <-  function(x){return(x[1]+(x[2]-x[1])/(1+exp((x[3]-cycles)/x[4])))}	 #define FPLM to search b with
    billies <- apply(bgrid,1,FPLM)							                         #generate possible outcomes
    bestb <- which.min(colSums(apply(billies,2,function(p){log(p)-log(fluo)})^2))  #best outcome is that which minimizes the difference between real and predict
    b <- bgrid[bestb,4]
	  sp<-c(Y0= Y0,a= a,X0= X0,b= b)
    #We now try to optimize starting values by perfroming a 4PLM fit (this is crucial, bad initial parameter estimates will screw with the 5PLM fit)
    #more so since the nls.lm routine has the tendency to keep soldiering on and return crazy parameter values instead of simply giving up
    mleko <- deriv(~Y0+(a-Y0)/(1+exp((X0-x)/b)),c("Y0","a","X0","b"),function(x,Y0,a,X0,b){})
	  modle<-try(nls(fluo~mleko(cycles,Y0,a,X0,b),start=sp,control=list(maxiter=100)),silent=TRUE)
	  #if the 4PLM fit fails, send message and use rough estimates
	  if(inherits(modle,"try-error")){
	       if(!isTRUE(silent)){message("Optimization of manual parameter values failed (4PLM fit)...")
                             message("using 'ad hoc' starting parameters.")}
	       }else{
	       sp<-coef(modle)}
   #final message (separator)
   if(!isTRUE(silent)){message("----")}
  }
  
##preliminary baseline subtraction (we use the 4plm fluo 5% fluo increase mark - 5 cycles as a first guess of where the groundphase ends) 
  #fit 4PLM (if it hasn't already been fit)
		stopg<-floor(sp[3]-log(0.95/0.05)*sp[4])-5	
      if(stopg<4){stopg<-8}#in case of failure, use minimal groundphase to prevent errors here. The usual error messaging will kick in downstream 
    sbase<-lm(fluo[4:stopg]~cycles[4:stopg])$coefficients
    sdata<-fluo-(sbase[2]*cycles+sbase[1])

##for debugging:   
  if(isTRUE(debugg)){
   x11(14,14);layout(matrix(c(1,4,2,2,3,3,5,5),nrow=4,ncol=2,byrow=F))
   plot(fluo[1:stopg]~cycles[1:stopg],bty="L",main="preliminary baseline");abline(sbase,col="red")}
    
##define model function (without slope for baseline)
	fRichL<-function(x,p){return(p[2]          + (p[1] - p[2]) / (1 + (2^(1/p[5])-1) * exp(p[4]*(x-p[3])))^p[5])}
##define Residual Sum of Squares function for nls.lm
  sRRs<-function(m){resi<-sdata-fRichL(cycles,m)}
##prep starting values for 5PLM (based on 4PLM parameters)    
  sp<-c(sp[1],sp[2],sp[3],(0.335-0.0195*sp[4]),1)
  names(sp)<-c("y0","fmax","xmid","b","g")
##fit 5PLM model to data
  sp<-try(nls.lm(par=sp,fn=sRRs),silent=TRUE)
  if(inherits(sp, "try-error")|any(is.nan(sp$par)==T)){
      sp<-c(NA,NA,NA,NA,NA)   #model parameters
      resec<-rep(NA,times=endd)  #model residuals
      fail<-201
      #Failure to fit a 5PLM may be a symptom of absence of amplification
      #if there is less than a 32 fold increase in fluorescence there is probably no amplification (less then 5 cycles of amplification)
      if((fluo[endd]/fluo[1])<32){fail<-202}                   
  }else{
  resec<-resid(sp)
  sp<-sp$par}
##add preliminary baseline to model parameters
  sp<-c(sp[1]+sbase[1],sbase[2],sbase[1]+sp[2],sp[3],sp[4],sp[5])  
  names(sp)<-c("y0","s","fmax","xmid","b","g")
##redefine model function (WITH slope for baseline)
  fRichL<-function(x,p){return(p[3] + p[2]*x + (p[1] - p[3]) / (1 + (2^(1/p[6])-1) * exp(p[5]*(x-p[4])))^p[6])}

##for debugging:   
  if(isTRUE(debugg)){plot(fluo~cycles,bty="L",main="5PLM fit");curve(fRichL(x,sp),0,endd,add=T)}

########################################
##          5PLM fit Check            ##
########################################

#We calculate the residual error per cycle + check the tolerance limit, if the [resec] value exceeds the limit we quit analysis
  if(all(is.na(resec))){#if all residuals are NA then the fit failed and we set resec to larger than its tolerance
    resec<-resec.tolerance*10
  }else{
    resec<-sum(abs(resec),na.rm=T)/endd}
##for debugging:   
  if(isTRUE(debugg)){message("resec value");print(resec);message("----")}

#if the [resec] value does not exceed its tolerance we may proceed  
if(!resec>=resec.tolerance){

#We warn the users in case of high values for [resec]
if(!isTRUE(silent)){
      if(resec>0.0675){message("Warning: exceptionally high residual error per cycle [resec value]!") ; message("data most likely very noisy");message("----")
      }else{if(resec>0.053){message("Warning: elevated residual error per cycle [resec value]!"); message("data may be noisy");message("----")}}}

#We calculate the necessary values from the 5PLM fit
  #Second derivative maxima aka the Cq value(s)
    sdm<-as.numeric((sp[5]*sp[4]+log(-(1/2)*(-1-3*sp[6]+sqrt(1+6*sp[6]+5*sp[6]^2))/(sp[6]^2*(2^(1/sp[6])-1))))/sp[5])
    sdm2<-as.numeric((sp[5]*sp[4]+log((1/2)*(1+3*sp[6]+sqrt(1+6*sp[6]+5*sp[6]^2))/(sp[6]^2*(2^(1/sp[6])-1))))/sp[5])
    #We use the y-axis position of the 5PLM sdm as threshold:
    ysdm<-fRichL(sdm,sp)-(sp[1]+sdm*sp[2]) #substract baseline!

  #calculate reference points of fluorescence accumulation
    fluospar<-function(p){return(floor(log(((-1/(p-1))^(1/sp[6])-1)/(2^(1/sp[6])-1))/sp[5]+sp[4]))}
    stada<-fluospar(0.075)           #cycle in which 7.5% of total fluorescence increase is reached
    stoda<-fluospar(0.85)             #cycle in which 85% of total fluorescence increase is reached
      #in case of data truncation the plateau normally follows directly after the truncation causing the 
      #95%FU mark to lie within the available cycle range, this may not always be the case
      #therefore it is useful to check if the 85% mark does not lie outside of the available cycle range
      if(stoda>endd){stoda<-endd}
    stp<-fluospar(0.95)               #cycle in which 95% of total fluorescence increase is reached
    startg<-4                         #start groundphase: skipping first three cycles to avoid high-noise data
    stopg <-stada-7                   #stop groundphase: skipping cycles that may be contaminated with amplicon fluorescence

##for debugging:   
  if(isTRUE(debugg)){message("Demarcation:");print(c("startg"=startg,"stopg"=stopg,"stada"=stada,"stoda"=stoda,"stp"=stp,"endd"=endd));message("----")
                              abline(v=c(startg,stopg,stada,stoda,stp),lty=2,col="snow4")}

#If any of these return NA that's a strong indication that the 
#data do not contain a valid PCR curve
if(any(is.na(c(stada,stoda,stp)))|any(c(stada,stoda,stopg,stp)<0)){fail<-203}

}else{fail<-204}#END of [resec] value check

#if the fit of the 5PLM failed OR the values calcuated from the fit are invalid
#there is no use in contiuing the procedure. Therefore we only proceed if fail equals zero
if(fail==0){
################################################################################
                ########################################
                ##         End of Pre-Analysis        ##
                ##              Start of              ##
                ##         Baseline subtraction       ##
                ########################################
################################################################################

##Baseline Function Repository
###################################
  #Restricted freedom model derivative
  mleko <- deriv(~-(a/b)^2-(a/b^2)*x-(2/b^2)*x^2,c("a","b"),function(x,a,b){})
  #BERT: quadratic model 
  Bert<-function(n,l){return( l[1]+n*l[2]+l[3]*n^2 )}   
  #HEMAN: shifts the baseline up or down a small amount and returns the sum of the E-residuals
  #       this will allow to find the basline shift that yields an expected value of the residuals = zero (irrespective of the size of the residuals!)
  heman<-function(cc,Rs=F){#set Rs to T to return residuals
                    Fmp<-sdata-cc
                    Lmp<-c(NA,log(abs(log(Fmp[-1]/Fmp[-endd]))))
                    #quick-fit of quadratic model and estimation of Emax
	                  pmp<-try(nls(Lmp~mleko(Fmp,a,b),start=lap,control=list(maxiter=100)),silent=TRUE)
                    if(inherits(pmp, "try-error")){Emp<-log(log(1))}else{pmp<-coef(pmp);Emp<--(pmp[1]/pmp[2])^2}
                    Rmp<-c(NA,log(abs(log(Fmp[-1]/Fmp[-endd]))))[startg:stopg]-Emp
                    #if Rs=T the residuals are all we need
                      if(isTRUE(Rs)){
                        return(Rmp)
                      }else{#we continue by exluding the outliers and calculating the residuals' sum
                        not<-sapply(boxplot.stats(Rmp)$out,function(q){return(which(Rmp==q))})
                        if(mode(not)=="list"){not<-unique(unlist(not))}#remove possible double values
                        if(length(not)==0){#there are no outliers
                          return(abs(sum(Rmp,na.rm=T)))
                        }else{            #there are outlies
                          return(abs(sum(Rmp[-not],na.rm=T)))
                        }}}
  #PEAKFINDER: the funcion that finds all local minima in a range of data
  findpeaks <- function(vec,bw=1,x.coo=c(1:length(vec))){  #where bw = is box width, setting the sensitivity of the search
                    ###set all vectors to null
                   	pos.x.max <- NULL ;	pos.y.max <- NULL ;	pos.x.min <- NULL ;	pos.y.min <- NULL
                    ###Start of for loop:    we walk down the vector with a window of size "bw"
                    for(i in 1:(length(vec)-1)){
                      #check if we have reached the end of the vector
                      if((i+1+bw)>length(vec)){sup.stop <- length(vec)}else{sup.stop <- i+1+bw}
                      #check if we are at beginning of the vector
                  		if((i-bw)<1){inf.stop <- 1}else{inf.stop <- i-bw}
                      #select window in two parts: values beyond i (superior), and values before i (inferior)
                  		subset.sup <- vec[(i+1):sup.stop]
                  		subset.inf <- vec[inf.stop:(i-1)]
                    ##############################################################
                      #are ALL trailing data smaller than i?
                  		is.max   <- sum(subset.inf > vec[i]) == 0
                  		#are ALL leading data smaller than i?
                   		is.nomin <- sum(subset.sup > vec[i]) == 0
                      #are ALL trailing data larger than i?
                   		no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
                      #are ALL leading data larger than i?
                  		no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
                    ##############################################################
                      #a maximum is found if  all data before and after i are smaller than i
                  		if(is.max & is.nomin){
                  			pos.x.max <- c(pos.x.max,x.coo[i])
                   			pos.y.max <- c(pos.y.max,vec[i])}
                      #a maximum is found if  all data before and after i are larger than i
                  		if(no.max & no.nomin){
                  			pos.x.min <- c(pos.x.min,x.coo[i])
                  			pos.y.min <- c(pos.y.min,vec[i])}}#end of for loop
                      ###Output
                       return(list("max.X"=pos.x.max,"max.Y"=pos.y.max,"min.X"=pos.x.min,"min.Y"=pos.y.min))} 
#Repository END####################

##Baseline Subtraction
###################################
##check for very "early" reactions (Cq<10, too few cycles for baseline fitting)
  if(base.line=="flat"|(stopg-startg)<1){
    if(!isTRUE(silent)&(stopg-startg)<1){message("Groundphase too short for linear model!");message("Using flat baseline instead...");message("----")}    
    sbase<-c(mean(fluo[startg:stopg]),0)
    }else{
##fit baseline (linear model: least squares, normal error).
    sbase<-lm(fluo[startg:stopg]~cycles[startg:stopg])$coefficients}
    names(sbase)<-c("intercept","slope")
    baseline<-sbase[2]*cycles+sbase[1]
    sdata<-fluo-baseline 

##Baseline optimisation
###################################
if(base.line=="optimized"){
#prep initial model estimate
   Fnsub<-sdata[stada:stoda]
   Lnsub<-c(NA,log(abs(log(sdata[-1]/sdata[-endd]))))[stada:stoda]
   lap<-coef(lm(Lnsub~Fnsub+I(Fnsub^2)))
   lap<-c(sqrt(abs(lap[1]))*sqrt(1/abs(lap[3])),sqrt(1/abs(lap[3])))
   names(lap)<-c("a","b")
#Let heman do the heavy lifting :                      
  temp<-seq(from=-0.005*sp[3],to=0.005*sp[3],length=200)   #we scan a region of 0.5% of the max fluo increase around zero with grid:200
  temp<-rbind(temp,sapply(temp,heman,Rs=F))                #we calculate the abs(sum) of the residuals for each value
  #detrending the values (helps to locate the relevant minima)
  temp<-rbind(temp,temp[2,]-(temp[2,1]-(((temp[2,200]-temp[2,1])/(temp[1,200]-temp[1,1]))*temp[1,1])+((temp[2,200]-temp[2,1])/(temp[1,200]-temp[1,1]))*temp[1,]))
#actual searching for minima:
  pilipili<-findpeaks(temp[3,])$min.X
  
##for debugging:   
    if(isTRUE(debugg)){plot(temp[3,]~temp[1,],bty="L",main="baseline shift");abline(v=temp[1,pilipili],lty=2,col="snow4")}

  #restrict possibilities to minima
  temp<-temp[,pilipili] 
  #now select lowest value
  if(is.null(dim(temp))){#there is only one minimum
    skeletor<-temp[2]
  }else{#multiple minima found
    skeletor<-temp[1,which.min(temp[2,])]}
    
##for debugging:   
    if(isTRUE(debugg)){abline(v=skeletor,lty=2,col="green3")
    plot(fluo[1:stopg]~cycles[1:stopg],bty="L",main="final baseline");abline(sbase,col="orange")}    

#update the baseline    
    if(!is.na(skeletor)){
    #optimisation succeeded
      sbase[1]<-sbase[1]-skeletor  #uptdate baseline paramters
        ##for debugging:
        if(isTRUE(debugg)){abline(sbase,col="green3");
                           message("baseline shift:");print(skeletor);message("----")
                           plot(c(NA,log(abs(log(sdata[-1]/sdata[-endd]))))~sdata,bty="L",main="ln E estimates",col="orange")}
      baseline<-baseline-skeletor #update baseline
      sdata<-sdata-skeletor       #update baseline subtracted fluorescence values
    }else{#If optimisation failed: we keep the initial estimate & notify user
      if(!isTRUE(silent)){message("Baseline optimisation failed!");message("Using initial estimate...");message("----")}}     
}#end of Baseline optimization if-module

##Now that we've obtained the baseline subtraced fluorescence values we may calculate the efficiency estimates                  
    Endata<-c(NA,(sdata[-1]/sdata[-endd]))
    #in Phase 2 the double log often causes NaNs which we cannot use in therefore we resort to
    #using abs() as a means of data imputation to counteract loss of points due to NaN from log(-)   
    Lndata<-c(NA,log(abs(log(sdata[-1]/sdata[-endd]))))
          
##for debugging:   
    if(isTRUE(debugg)){points(Lndata~sdata,pch=2,col="green3")
                      abline(v=sdata[c(startg,stopg,stada,stoda,stp)],lty=2,col="snow4")}
    
################################################################################     
                ########################################
                ##         End of Pre-Analysis        ##
                ##             Start of               ##
                ##     General Function Repository    ##
                ########################################
################################################################################     

  #Define the restricted model
    Billy<-function(n,l){return(-(l[1]/l[2])^2-(l[1]/l[2]^2)*n-(2/l[2]^2)*n^2)}
  #Define bilinear model
    Bob<-function(n,b){return(b[6]+b[5]*log(exp((b[1]*(n-b[4])^2+b[2]*(n-b[4]))/b[5])+exp(b[3]*(n-b[4])/b[5])))}                          
  #Function to calculate the initial fluorescence (a*i0)  
    iniFlu<-function(kp,kb){
         ##In the baseline subtracted data we replace the groundphase by an estimate of the baseline subtracted fluo based on the Emax we found
         ##this prevents the noise in the baseline from affecting the ai0 estimate
         #first copy baseline subtracted data
         sdataplus <- sdata
         for(i in stada:1){
         sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],kp)))} #Bert will do for both phasees since only the very beginning of the first phase is involved
         #define weights for ai0 determination process
         weight<-1.0001-sdataplus/sdataplus[endd]
         ##The upper and lower limit for the initial target copies * amplicon Fluorescence are  1e+10 and 1e-10 respectivly
         opti<-optimize(mini_one,interval=c(1e-10,1e+10),w=weight,b=kb,splus=sdataplus)
         #we multiply the outcome with the scaling factor for the final result
         result<-opti$minimum*1e-10
         return(result)}   
  #Function to calculate the BCa interval from the bootstrap estimates (bot), jackknife estimates (jac) & the original sample estimate (ori)
    BCa<-function(ori,bot,jac,etna="parameter"){
         #check for NA values in bootstrap and jackknive estimates. Function is robust for NAs but we should at least warn user.
         if(any(is.na(bot))&!isTRUE(silent)){message(paste(etna,"bootstrap estimates contain NA value(s)",sep=" "))}
         if(any(is.na(bot))&!isTRUE(silent)){message(paste(etna,"jackknive estimates contain NA value(s)",sep=" "))}
         #in order to prevent z = - inf we check if the denominator is zero before calculatig z
         if(length(which(bot<=ori))==0){z<-qnorm(0.1/(n+1))}else{z<- qnorm(length(which(bot<=ori))/(n+1))}                                                  
         a<- (sum((jac-mean(jac,na.rm=T))^3,na.rm=T))/(6*sum((jac-mean(jac,na.rm=T))^2,na.rm=T)^(3/2))   
         a1<- pnorm(z+(z-qnorm(1-p/2))/(1-a*(z-qnorm(1-p/2))))                                           
         a2<- pnorm(z+(z+qnorm(1-p/2))/(1-a*(z+qnorm(1-p/2))))                                          
         upp <- sort(bot)[ceiling(n*a1)] #using 'ceiling' rather than 'floor' or 'round' to prevent asking for index 'zero'
         low <- sort(bot)[ceiling(n*a2)] 
         BCout<-c(upp,low)
         names(BCout)<-c("upper","lower")
         return(BCout)}
  #Function to generate the (double log) cycle efficiencies from a set of parameters
    Lngen<-function(kp,kb=NA,ph=1){
        sdataplus <- sdata
        #now replace groundphase cycles with kinetic estimates
        for(i in stada:1){
        sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],kp)))}
        #calculate individual cycle eff estimates
        if(ph==1){LnSet<-Bert(sdataplus,kp)}
        if(ph==2){LnSet<-Bob(sdataplus,kb)}
        return(LnSet)}
  #Function to find the Cq value given a threshold (using spline function): 
    get.Ct<-function(fluo,threshold,RCt){#Rct is a rough Ct estimate (the spline will search for a Ct in a 3 cycle radius)
        cycles<-c(1:length(fluo))
        smother<-try(splinefun(cycles,fluo,method="fmm"),silent=T)
        if(inherits(smother,"try-error")){ #For currently unknown reasons the construction of the splinefunction can go haywire (Error in splinefun(cycles, fluo, method = "fmm") : zero non-NA points)
          Ct<-NA                           #In order not to interrupt the entire workflow we return NA + a warning message
          if(!isTRUE(silent)){message("NA value during Ct calculation")}
        }else{
          Ct<-optimize(function(x){return(abs(threshold-smother(x)))},interval=c(RCt-3,RCt+3))$minimum}  #limit search to a small window to avoid local minima
        return(Ct)}

################################################################################     
                ########################################
                ##          End of Repository         ##
                ##    Start of actual Model Fitting   ##
                ########################################
################################################################################     

########################################
##     First Phase Calculations       ##
########################################

#Irrespective of the number of phases in the data, we alays need the first phase of decline
#Therefore, we fit a linear model to first phase

##select data cycles & respective fluorescence values
   Fnsub<-sdata[stada:stoda]
   Lnsub<-Lndata[stada:stoda]
##for debugging:   
  if(isTRUE(debugg)){
   x11(14,14);layout(matrix(c(1,1,2,2,3,3,4,4),nrow=4,ncol=2,byrow=F))
   plot(Lndata~sdata,bty="L",main="first phase fitting")}
   
##Start of the fitting story
    lpM<-try(lm(Lnsub~Fnsub+I(Fnsub^2)),silent=T)
  ##If the fit didn't work we take not and stop here
    if(inherits(lpM, "try-error")){
          fail<-append(fail,301)
  ##If the fit did work, we extract the model parameters
    }else{
  ##First we extract the model parameters
  lp<-coef(lpM) 
  #now we transform them into estimates for the restricted freedom model
  lap<-c(sqrt(abs(lp[1]))*sqrt(1/abs(lp[3])),sqrt(1/abs(lp[3])))
  names(lap)<-c("a","b")
    ##for debugging:   
    if(isTRUE(debugg)){
      curve(Bert(x,lp),0,sdata[stoda+1],add=T,lty=2,col="red")
      curve(Billy(x,lap),0,sdata[stoda+1],add=T,lty=2,col="orange")}    
  ##actual fit of the model:    
	  modle<-try(nls(Lnsub~mleko(Fnsub,a,b),start=lap,control=list(maxiter=100)),silent=TRUE) 	  
    ##check outcome
    if(inherits(modle,"try-error")){
	       if(!isTRUE(silent)){message("restricted freedom model failed");message("using all degrees of freedom....")}
    ##If the fit did work, we continue by applying the KALMAN filter
    }else{        
	  #First we extract and convert the parameters from the model
         lap<-coef(modle)
	       #convert to 3 parameter model
         lp<-c(-(lap[1]/lap[2])^2,-(lap[1]/lap[2]^2),-(2/lap[2]^2))
        ##for debugging:   
        if(isTRUE(debugg)){curve(Bert(x,lp),0,sdata[stoda+1],add=T,lty=1,col="green3")
                           plot(Lndata[1:stoda]~sdata[1:stoda],bty="L",main="Kalman Filter")
                           message("parameters before Kalman Filter:");print(lp)}

##If desired we continue by applying the KALMAN filter  
if(isTRUE(kalman)){         
  ##We calculate fluorescence increases between each cycle, we'll need these since the Kalman filter uses En-1 = En - Vn * deltaF + C * deltaF^2
    #rather than En-1 = E0 + B * Fn-1 + C * (Fn-1)^2
    deltaFs<-c(NA,abs(sdata[1:(stada-1)]-sdata[2:stada]))       
    
  ##Calculate the initial states for the filter to start from
    #current state
    currstate<-matrix(c(Bert(sdata[stada],lp),lp[2]+2*lp[3]*sdata[stada]),nrow=2,ncol=1)
    states<-currstate
    #we can calculate sd of the acceleration a using the residuals from the regression  (matrix form: sigma^2 =  MSE * (X'X)^-1, see page 207 of book)
      xix<-matrix(c(rep(1,times=length(Fnsub)),Fnsub,Fnsub^2),nrow=length(Fnsub),ncol=3,byrow=F)   
      MSE<-sum((Lnsub-Bert(Fnsub,lp))^2)/(length(Fnsub)-3)   #we estimted three parameters so we loose 3 df
      sig<- 2 * (MSE*diag(solve(t(xix)%*%xix)))[3]^(1/2)     # * 2 since our parameter c equals 2 a (acceleration) in the physical model
    #variance-covariance matrix of the estimates
    estcov<-matrix(c((sig*deltaFs[stada]^2)^2,sig^2*2*deltaFs[stada]^3,sig^2*2*deltaFs[stada]^3,(sig*2*deltaFs[stada])^2),nrow=2,ncol=2,byrow=T)

  ##we also calculate a measure of variance (sd) for each of the Ln(Ln(E)) estimates using a bootstrap approach
    #we generate a first sdataplus
         sdataplus <- sdata
         for(i in stada:1){
         sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],lp)))} #Bert will do for both phasees since only the very beginning of the first phase is involved
    #estimate noise using the groundphase
        sdb<-sd(resid(lm(fluo[startg:stopg]~cycles[startg:stopg]))) #since we shifted the baseline we can no longer us it to generate the residuals with E=0
    #we now make the sdataplus based bootstraps
    #it doesn't matter that the model is not 100% correct yet, it's mainly the influence of the noise on the estimates we're interested in, not the estimates themselves
      strps<-matrix(rep(sdataplus,times=1000)      ,nrow=1000,ncol=endd,byrow=T)  +  #fitted values
             matrix(rnorm(1000*endd,mean=0,sd=sdb) ,nrow=1000,ncol=endd)             #residuals
      strps<-cbind(rep(NA,times=1000),strps[,-1]/strps[,-endd])
      strps<-apply(log(log(strps)),2,sd,na.rm=T) #thus we have the measurement uncerntainty of the log(log(E))s

  ##Actual Filter
  ################
  for(i in stada:2){
  ##coefficient matrices
    #state prediction matrices
    AA<-matrix(c(1,-deltaFs[i],0,1),nrow=2,ncol=2,byrow=T)
    BB<-matrix(c(deltaFs[i]^2,-2*deltaFs[i]),nrow=2,ncol=1)
    #measurement prediction matrix
    CC<-matrix(c(1,0),nrow=1,ncol=2)
  ##Covariance matrices
    #state covariance matrix
    EE<-matrix(c((sig*deltaFs[i]^2)^2,sig^2*2*deltaFs[i]^3,sig^2*2*deltaFs[i]^3,(sig*2*deltaFs[i])^2),nrow=2,ncol=2,byrow=T)
    #measurement coveriance matrix (look up in bootstraps vector)
    EZ<-strps[i]   
  ##The run
    #state prediction
    statepred<-AA%*%currstate+BB%*%lp[3]
    #measurement prediction
    measupred<-CC%*%statepred
    #estimate covariance
    estcov<-AA%*%estcov%*%t(AA)+EE
    #Kalman gain
    K<-estcov%*%t(CC)%*%solve((CC%*%estcov%*%t(CC)+EZ))
    #Final state prediction
    statepred<-statepred+K%*%(Lndata[i-1]-measupred)
    states<-cbind(states,statepred)
    #Final covariance estimate
    estcov<-(diag(2)-K%*%CC)%*%estcov
  ##refit model
  ##refit model
    Fnsub<-sdata[(i-1):stoda]
    Lnsub<-c(rev(states[1,]),Lndata[(stada+1):stoda])
	  modle<-nls(Lnsub~mleko(Fnsub,a,b),start=lap,control=list(maxiter=100))
    #extract paramters
    lpK<-coef(modle)
    lpK<-c(-(lpK[1]/lpK[2])^2,-(lpK[1]/lpK[2]^2),-(2/lpK[2]^2))
      ##for debugging:   
      if(isTRUE(debugg)){points(Lnsub~Fnsub,col="orange");curve(Bert(x,lpK),0,sdata[stoda+1],add=T,lty=2,col="orange")}
    #recalculate sig
    xix<-matrix(c(rep(1,times=length(Fnsub)),Fnsub,Fnsub^2),nrow=length(Fnsub),ncol=3,byrow=F)
    MSE<-sum((Lnsub-Bert(Fnsub,lpK))^2)/(length(Fnsub)-3)  
    sig<- 2 * (MSE*diag(solve(t(xix)%*%xix)))[3]^(1/2)    
  ##advance current state using updated model
    currstate<-matrix(c(Bert(sdata[i-1],lpK),lpK[2]+2*lpK[3]*sdata[i-1]),nrow=2,ncol=1)
    }#END of the filter for-loop
  #we check if the last Kalman model is valid and, if so, update the model paramters
  if(!any(is.na(lpK))){
    lp<-lpK
      ##for debugging:   
      if(isTRUE(debugg)){curve(Bert(x,lp),0,sdata[stoda+1],add=T,lty=1,col="green3")
                        message("parameters after Kalman Filter:");print(lp)}
  }else{
    #if the filter failed we inform the user
    if(!isTRUE(silent)){message("Kalman filter failed, using initial model estimate...")}}
  }}}#END of the fitting story we either have Fail!=0 or lp exists

################################################################################
# we should only continue of we have parameters (lp exists, i.e. the first phase fit succeeded)!
if(exists("lp")){
################################################################################
  ##                   ##
  ## Bootstrap Block I ##
  ##                   ##
#########################
if(isTRUE(bootstrap)){
    #we generate the Residual Sum of Squares (residuals are Ln, the horizontal distances are no more )
    RSS1<-c(sum((Lndata[1:(stada-1)]-Bert(sdata[1:(stada-1)],lp))^2,na.rm=T),sum((Lndata[stada:stoda]-Bert(sdata[stada:stoda],lp))^2))
    if(isTRUE(kalman)){temp<-phaseboot(Fn=sdata,Ln=Lndata,dmc=c(startg,stopg,stada,stoda,endd),KAL=TRUE,Evar=strps,RSS=RSS1,lice=debugg,block=1,ose=lp,silent=T)
                 }else{temp<-phaseboot(Fn=sdata,Ln=Lndata,dmc=c(startg,stopg,stada,stoda,endd)                    ,RSS=RSS1,lice=debugg,block=1,ose=lp,silent=T)}
    #split
    Blps<-temp$blp 
    Jlps<-temp$jlp }

################################################################################
## IF no truncation can be detected we continue with the second phase fit & construction of the bilinear model
if(stp<(endd-4)){
################################################################################

## First: Fit the second phase of decline             
################################################################################     

##make data subset for fitting the second phase of effiency decline
    #start of second phase of efficiency decline == where fluorescence reaches 95% (as calculated from the 5PLM, see above)
    #select data cycles & respective fluorescence values
    Fnsub<-sdata[stp:endd]
    #in Phase 2 the double log often causes NaNs which we cannot use in therefore we resort to
    #using abs() as a means of data imputation to counteract loss of points due to NaN from log(-)
    Lnsub<-Lndata[stp:endd]  
    #perform model fit (linear, least squares, normal error) on the inverse relation 
    #Fn vs E, rather then E vs Fn
    #Fn is measured and holds the error, so the residual distance to be minized is the distance to Fn
    lpM2<-lm(Fnsub~Lnsub)
    lp2 <-coef(lpM2)          #extract parameters
    RSS2<-sum(resid(lpM2)^2)  #generate residual sum of squares (bootstrap purposes)
      #the probability exists that the linear model wil yield a positive slope
      #which causes a contradiction in the bilinear model
      #in those cases we force a negative slope through a custom fitting approach
      if(lp2[2]>=0){
      B.negative<-function(zz){return(sum((Fnsub-(zz[1]-zz[2]^2*Lnsub))^2,na.rm=TRUE))}
      lpM2<-optim(par=c(Fnsub[1],0),fn=B.negative)
      lp2 <-lpM2$par
      RSS2<-lpM2$value  #extract RRS2 for bootstrap purposes (in this case value of 'fn' at 'par', see optim) 
      lp2[2]<--lp2[2]^2}#update model parameter to correct value
    #convert the parameters for so the plot remains E vs Fn
    lp2<-c((-lp2[1]/lp2[2]),(1/lp2[2]))  
    ##for debugging:   
    if(isTRUE(debugg)){
       plot(Lndata~sdata,bty="L",main="second phase fitting")
       curve(Bert(x,lp),0,sdata[endd],add=T,lty=2,col="orange")
       abline(lp2,col="orange",lty=2)}
        
## Second: construct the bilinear model using the linear fits from both phases             
################################################################################     

  ##We first calculate the Emax (intercept) 
    Lmax<-lp[1]
    Emax<-exp(exp(lp[1]))
    if(isTRUE(bootstrap)){                                #for future use: when using reference Efficiency in stead of max efficiency we will need a highly similar code
    BEmax<-exp(exp(Blps[1,]))                             #BEmax<-exp(exp(apply(Bbps,2,Bob,n=0)))
    JEmax<-exp(exp(Jlps[1,]))                             #JEmax<-exp(exp(apply(Jbps,2,Bob,n=0)))
    BLmax<-Blps[1,]                                       #BLmax<-apply(Bbps,2,Bob,n=0)
    JLmax<-Jlps[1,]}                                      #JLmax<-apply(Jbps,2,Bob,n=0)}

  #then we estimate where both phases cross 
    crx <- (1/2)*(lp2[2]-lp[2]+sqrt(lp2[2]^2-2*lp2[2]*lp[2]+lp[2]^2+4*lp[3]*lp2[1]-4*lp[3]*lp[1]))/lp[3] 
    
  ##finally, we gather all intial parameters for the bilinear model
    slo1a<-lp[3]
    slo1b<-lp[2]+2*crx*lp[3]    #see remark:
        #for simple linear components of the bilinear model, the slope is not affected by the (x-crx) term (displacement of the function along the x-axis).
        #Only the intercept changes when the function is displaced by 'crx', and the intercept is not used in the bilinear model.
        #For quadratic equations the second slope (b) and intercept (int) change by displacement along the x-axis: int + b*(x+crx) + a*(x+crx)^2 <=> (int+b*crx+a*crx^2) + x * (b+2*crx*a) + a*x^2
        #thus displacement of the curve by 'crx' changes the 'b' slope, hence we adjust its value.
    slo2<-lp2[2]  
    eta<--4
    chi<-lp[1]-eta*log(exp((slo1a*crx^2-slo1b*crx)/eta)+exp(-slo2*crx/eta))
    bp<-c(slo1a,slo1b,slo2,crx,eta,chi)

##Only eta should be optimized. Shifting crx will shift the whole of the second phase horizontally left or right
    #redifine action zone
    Fnsub<-sdata[(stoda-2):endd]
    Lnsub<-Lndata[(stoda-2):endd]      

    bobslee<-function(d){ #use E values to make it more linear? and the RSS larger
      return(
      sum((exp(exp(Lnsub))-
         exp(exp((lp[1]-d*log(exp((bp[1]*bp[4]^2-bp[2]*bp[4])/d)+exp(-bp[3]*bp[4]/d))
         +d*log(exp((bp[1]*(Fnsub-bp[4])^2+bp[2]*(Fnsub-bp[4]))/d)+exp(bp[3]*(Fnsub-bp[4])/d)))))      
          )^2,na.rm=TRUE)
          )}
    crash<-try(optim(fn=bobslee,par=bp[5],method="L-BFGS-B",lower=-10,upper=-2.5),silent=TRUE)
    #this optimalisation may crash, in that case we keep the original values and take notice of the error
    if(inherits(crash, "try-error")){
    fail<-append(fail,307)
    }else{
    bp[5]<-crash$par
   #update Chi
    bp[6]<-lp[1]-bp[5]*log(exp((bp[1]*bp[4]^2-bp[2]*bp[4])/bp[5])+exp(-bp[3]*bp[4]/bp[5]))}
    names(bp)<-c("a2","a1","a3","ipt","eta","chi")   
    ##for debugging:   
    if(isTRUE(debugg)){curve(Bob(x,lp),0,sdata[endd],add=T,lty=2,col="green3")}

################################################################################
  ##                    ##
  ## Bootstrap Block II ##
  ##                    ##
##########################
if(isTRUE(bootstrap)){
    Fnsub<-sdata[stp:endd]
    Lnsub<-Lndata[stp:endd]
    #we now hand all material to the bootstrap function. Note that to be able to calculate the crx-point
    #the function needs to know the lm-paramters or the previous phase, thus we also input Bpr1 & Jpr1.
    #Finally, to complete the full set of trilinear jackknives we also need the original sample parameters, hence the input of  lp1 and lp2
    temp<-phaseboot(Fn=Fnsub,Ln=Lnsub,RSS=RSS2,Bpar1=Blps,Jpar1=Jlps,ose=c(lp,lp2,bp),block=2,silent=FALSE)
    #split
    Bbps<-temp$bbp  
    Jbps<-temp$jbp}

## Third: define the fucntion to rebuild the fluorescence readings using the calculated cycle efficiencies  (TWO phase version)
################################################################################     

##Here comes the Residual sum of squares function
##it takes the 'fluorescence due to initial target copies' as input value since that is the value we want to optimize
    mini_one<-function(m,w,b,splus){return(
        #we'll calculate the sum of squared residuals for the holistic model, so we'll start with the sum command 
				sum(w*                                                                    
        #we substract the data with our own Fn calculations to yield the model residuals
				(log(fluo) -                                                                  
        #begin of cycle Fluorescence calculation
        log(                                                                     
          ####baseline (duh)
					baseline +                                                                   
          ###Fluorescence due to initial target copy input (alpha*i_0)
           ## the '1e-10' is a scaling factor, which helps the algorithm 
           ## to reach very small values
          m *  1e-10  *
          ###The cumumlative product of reaction efficiencies: this will yield a value for each cycle.
          ###Due to all of the above this value is converted in a fluorescence estimate
          ###which in turn yields the residuals after subtraction from the actual data
          cumprod(exp(exp(Bob(splus,b))))                                                            
        #end of cycle Fluorescence calculation
        )                                                        
				#now we just have to sqaure the residuals and sum them up
        )^2)
    ##that's all folks
    )}

########################################
##        Estimation of a*i0          ##
########################################

    ai0<-iniFlu(kp=lp,kb=bp)
    Lns<-Lngen(kp=lp,kb=bp,ph=2)
    Ens<-(exp(exp(Lns)))
    Fns<-baseline + ai0 * cumprod(Ens)                                                        

    if(isTRUE(debugg)){
         sdataplus <- sdata
         for(i in stada:1){
         sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],lp)))} #Bert will do for both phasees since only the very beginning of the first phase is involved
         #define weights for ai0 determination process
         weight<-1.0001-sdataplus/sdataplus[endd]
         alfis<-vector()
         for(i in 1:endd){
         ##The upper and lower limit for the initial target copies * amplicon Fluorescence are  1e+10 and 1e-10 respectivly
         opti<-optimize(function(m){return(sum(weight[1:i]*(log(fluo[1:i])-log(baseline[1:i]+m*1e-10*cumprod(exp(exp(Bob(sdataplus[1:i],bp))))))^2))},interval=c(1e-10,1e+10))
         #we multiply the outcome with the scaling factor for the final result
         alfis<-append(alfis,opti$minimum*1e-10)}
         plot(alfis,bty="L",main="ai0 per cycle",log="y")}
    
    if(isTRUE(bootstrap)){    
    Bi0s<-sapply(c(1:n),function(a){return(iniFlu(Blps[,a],Bbps[,a]))})         
    #for the jackknives things are somewhat less straightforward: to cover the sequential omision of all datapoints
    #for the construcion of the  biliniear model we should also make the first phase jackknives into bilinear models
    #(using the ose 2nd phase), to that we attach the trilinear models form the second phase jackknives (using the ose 1st phase)
    temp<-cbind(Jlps,matrix(rep(lp,times=(length(Jbps[1,])-length(Jlps[1,]))),nrow=3,ncol=(length(Jbps[1,])-length(Jlps[1,])),byrow=F))     
    Ji0s<-sapply(c(1:length(Jbps[1,])),function(a){return(iniFlu(temp[,a],Jbps[,a]))})                                                
    #next we generate the log(log(E))-chains for the bootstraps and jackknives
    Blns<-sapply(c(1:n),function(a){return(Lngen(Blps[,a],Bbps[,a],ph=2))})                                                           
    Jlns<-sapply(c(1:length(Jbps[1,])),function(a){return(Lngen(temp[,a],Jbps[,a],ph=2))})                                            
    }	   

################################################################################
#In case no second phase could be detected:
}else{ #part of the data truncation structure 
################################################################################
  #Inform user that data are truncated and no plateau was detected
    if(!isTRUE(silent)){message("Failed to detect second phase");message("Analysis limited to first phase of efficiency decline...");message("----")}
    
  #we calculate the maximal efficiency values (in double log and standard)
    Lmax<-lp[1]
    Emax<-exp(exp(lp[1]))
    if(isTRUE(bootstrap)){                                #for future use: when using reference Efficiency in stead of max efficiency we will need a highly similar code
    BEmax<-exp(exp(Blps[1,]))                             #BEmax<-exp(exp(apply(Bbps,2,Bob,n=0)))
    JEmax<-exp(exp(Jlps[1,]))                             #JEmax<-exp(exp(apply(Jbps,2,Bob,n=0)))
    BLmax<-Blps[1,]                                       #BLmax<-apply(Bbps,2,Bob,n=0)
    JLmax<-Jlps[1,]}                                      #JLmax<-apply(Jbps,2,Bob,n=0)}
  #we also change the plotting output, so something useful can be produced
    plotto<-2

## define the fucntion to rebuild the fluorescence readings using the calculated cycle efficiencies  (ONE phase version)
################################################################################     

  ##the incomplete data case requires a dedicated Residual sum of squares function
  ##it takes the 'fluorescence due to initial target copies' as input value since that is the value we want to optimize.
  ##Also: it does not use the data beyond the 85% fluorescence mark since the since validity of the efficiency model
  ##beyond that point is not guaranteed
    mini_one<-function(m,w,b,splus){return(
        #we'll calculate the sum of squared residuals for the holistic model, so we'll start with the sum command
				sum(w[1:stoda]*
        #we substract the data with our own Fn calculations to yield the model residuals
				(log(fluo[1:stoda]) -
        #begin of cycle Fluorescence calculation
        log(
          ####baseline (duh)
					baseline[1:stoda] +
          ###Fluorescence due to initial target copy input (alpha*i_0)
           ## the '1e-10' is a scaling factor, which helps the algorithm
           ## to reach very small values
          m *  1e-10  *
          ###The cumumlative product of reaction efficiencies: this will yield a value for each cycle.
           ##Due to all of the above this value is converted in a fluorescence estimate
           ##which in turn yields the residuals after subtraction from the actual data
          cumprod(exp(exp(Bert(splus[1:stoda],b))))
        #end of cycle Fluorescence calculation
        )
				#now we just have to sqaure the residuals and sum them up
        )^2)
    ##that's all folks
    )}
 
########################################
##        Estimation of a*i0          ##
########################################
	    
    ai0<-iniFlu(kp=lp,kb=lp)
    Lns<-Lngen(kp=lp,ph=1)
    Ens<-(exp(exp(Lns)))
    Fns<-baseline + ai0 * cumprod(Ens)                                                        

    if(isTRUE(debugg)){
         sdataplus <- sdata
         for(i in stada:1){
         sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],lp)))} #Bert will do for both phasees since only the very beginning of the first phase is involved
         #define weights for ai0 determination process
         weight<-1.0001-sdataplus/sdataplus[endd]
         alfis<-vector()
         for(i in 1:stoda){
         ##The upper and lower limit for the initial target copies * amplicon Fluorescence are  1e+10 and 1e-10 respectivly
         opti<-optimize(function(m){return(sum(weight[1:i]*(log(fluo[1:i])-log(baseline[1:i]+m*1e-10*cumprod(exp(exp(Bert(sdataplus[1:i],lp))))))^2))},interval=c(1e-10,1e+10))
         #we multiply the outcome with the scaling factor for the final result
         alfis<-append(alfis,opti$minimum*1e-10)}
         plot(alfis,bty="L",main="ai0 per cycle",log="y")}

    if(isTRUE(bootstrap)){                                                                     
    Bi0s<-sapply(c(1:n),function(a){return(iniFlu(Blps[,a],Blps[,a]))})                         
    Ji0s<-sapply(c(1:length(Jlps[1,])),function(a){return(iniFlu(Jlps[,a],Jlps[,a]))})          
    #next we generate the log(log(E))-chains for the bootstraps and jackknives
    Blns<-sapply(c(1:n),function(a){return(Lngen(Blps[,a],ph=1))})                              
    Jlns<-sapply(c(1:length(Jlps[1,])),function(a){return(Lngen(Jlps[,a],ph=1))})               
    }  

   
################################################################################    
}#END of the data truncation structure (from here on, all commands are valid irrespective of the number of phases in the data)
################################################################################

Ct<-get.Ct(sdata,threshold=ysdm,RCt=sdm)            #Origincal Ct

    ##for debugging:   
    if(isTRUE(debugg)){
       x11()
       plot(fluo~cycles,bty="L",main="ai0 fitting")
       points(Fns,col="green3",pty=16,cex=0.8)
       lines(Fns,col="green3")}



########################## 
  ##                    ##
  ## Bootstrap Block III##
  ##                    ##
########################## 
##This block calculates the various confidence intervals

if(isTRUE(bootstrap)){ 
##confidence interval for the ai0 estimate
  CIai0<-BCa(ai0,bot=Bi0s,jac=Ji0s,etna="ai0")
##confidence interval for the Emax estimate
  CIEmax<-BCa(Emax,bot=BEmax,jac=JEmax,etna="Emax")
##confidence interval for the parameter estimates
  #note that we will have to use different command depending on data completeness (bilinear / trilinear)
  if(exists("Bbps")){
    #if Bbps exists a Bilinear model has been fitted so we calculate CI for both the first phase and the bilinear model
    CIlip<-sapply(c(1:3),function(a){return(BCa(lp[a],bot=Blps[a,],jac=Jlps[a,],etna="quadratic parameters"))})
    CIbip<-sapply(c(1:6),function(a){return(BCa(bp[a],bot=Bbps[a,],jac=Jbps[a,],etna="Bilinear parameters"))})
  }else{
    #if bps does not exist we calculate CI for the single phase parameters only
    CIlip<-sapply(c(1:3),function(a){return(BCa(lp[a],bot=Blps[a,],jac=Jlps[a,],etna="quadratic parameters"))})
    #and we return a matrix of NAs in order to have a uniform output format
    CIbip<-matrix(rep(NA,times=12),nrow=2,ncol=6)}
##The individual cycle fluorescence & efficiency confidence intervals
  #we treat every cycle/fluorescence as a separate bootstrap, calculating as many confidence intervals as there are cycles/fluorescence
  #For the fluorescence values we re-add the baseline
  BFns<-sapply(c(1:n)               ,function(x){return(baseline+cumprod(exp(exp(Blns[,x])))*ai0)})
  JFns<-sapply(c(1:length(Jlns[1,])),function(x){return(baseline+cumprod(exp(exp(Jlns[,x])))*ai0)}) 
  #now that we have the Fluo values we can also calculate the Cts
  BCts<-apply(BFns,2,get.Ct,threshold=ysdm,RCt=sdm) #Bootstrap Cts
  JCts<-apply(JFns,2,get.Ct,threshold=ysdm,RCt=sdm) #Jackknive Cts
##confidence intervals for each cycle's Efficiency estimate & Fluorescence estimate
  CIFns<-sapply(cycles,function(a){return(BCa(Fns[a],bot=BFns[a,],jac=JFns[a,],etna="Fluorescence"))})
  CIEns<-sapply(cycles,function(a){return(BCa(Ens[a],bot=exp(exp(Blns[a,])),jac=exp(exp(Jlns[a,])),etna="Cycle Efficiency"))})
##confidence interval for each threshold cycle
  CICt<-BCa(Ct,bot=BCts,jac=JCts,etna="Ct")  
  }

################################################################################    
                ########################################
                ##         End of Calculations        ##
                ##       Start of Output Section      ##
                ########################################
################################################################################     

########################################
##         Plotting Section           ##
########################################

#Should there be plots?
if(isTRUE(plots)){

  ##Plot A: fluorescence plots
  #######################################  
   
   #open graphical device
    x11(width=14,height=7)
    par(mfcol=c(1,2))

   #depending on analysis print different plot:
   
  ##plot if plateau is reached:
  if(plotto==1){

    #optional plot 1: fitted model on normal scale
    plot(fluo~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",bty="L")
    points(Fns,col="green3",pty=16,cex=0.8)
    lines(Fns,col="green3")
    if(isTRUE(bootstrap)){
    #confidence interval
    lines(cycles,CIFns[1,],col="gold",lty=2);lines(cycles,CIFns[2,],col="gold",lty=2)}
    #mark the data partitions
    lines(c(startg,startg),c(sbase[1],(sp[3]/15)),lty=2)            #start baseline calculation window
    lines(c(stopg,stopg),c(sbase[1],(sp[3]/15)),lty=2)  #stop baseline calculation window
    lines(c(stada,stada),c((fluo[stada]-(sp[3]/20)),(fluo[stada]+(sp[3]/20))),lty=3,col="blue") #start first decline phase
    lines(c(stoda,stoda),c((fluo[stoda]-(sp[3]/20)),(fluo[stoda]+(sp[3]/20))),lty=3,col="blue") #stop first decline phase
    lines(c(stp,stp),c((fluo[stp]-(sp[3]/20)),(fluo[stp]+(sp[3]/20))),lty=3,col="red") #start second decline phase
    mtext("A",side=3,cex=1.5,adj=1)
    
    #optional plot 2: fitted model on log scale
    plot(sdata~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",log="y",bty="L")
    points((Fns-baseline),col="green3",pty=16,cex=0.8)
    lines((Fns-baseline),col="green3")
    if(isTRUE(bootstrap)){
      #confidence interval
      lines(cycles,(CIFns[1,]-baseline),col="gold",lty=2);lines(cycles,(CIFns[2,]-baseline),col="gold",lty=2)
      #write some text in margin as legenda:
      mtext("original sample estimate",col="green3",cex=0.8,side=1,line=2,adj=-0.15)
      mtext("confidence limits",col="gold",cex=0.8,side=1,line=3,adj=-0.13)}   
    abline(h=ysdm) #plot threshold
    if(isTRUE(bootstrap)){
      #indicate Ct confidence interval
      lines(c(Ct,Ct),c(-1000,fluo[floor(Ct)]),lty=3,col="green3")
      lines(c(CICt[1],CICt[1]),c(-1000,CIFns[1,floor(CICt[1])]),lty=3,col="gold")
      lines(c(CICt[2],CICt[2]),c(-1000,CIFns[2,floor(CICt[2])]),lty=3,col="gold")}
    mtext("B",side=3,cex=1.5,adj=1)
    }
    
  ##plots if plateau is NOT reached    
  if(plotto==2){

    #optional plot 1: fitted model on normal scale
    plot(fluo~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",bty="L")
    points(Fns,col="green3",pty=16,cex=0.8)
    lines(Fns,col="green3")
    if(isTRUE(bootstrap)){
    #confidence interval
    lines(cycles,CIFns[1,],col="gold",lty=2);lines(cycles,CIFns[2,],col="gold",lty=2)}
    #mark the data partitions
    lines(c(startg,startg),c(sbase[1],(sp[3]/15)),lty=2)            #start baseline calculation window
    lines(c(stopg,stopg),c(sbase[1],(sp[3]/15)),lty=2)  #stop baseline calculation window
    lines(c(stada,stada),c((fluo[stada]-(sp[3]/20)),(fluo[stada]+(sp[3]/20))),lty=3,col="blue") #start first decline phase
    ##mark the part of the data that was not used for fitting with red
    rect((stoda+0.1),-500,65,(max(fluo)+1000),col="coral",border=NA)  
    abline(v=stoda,lty=2)
    points(fluo~cycles)
    mtext("A",side=3,cex=1.5,adj=1)

    #optional plot 2: fitted model on log scale
    plot(sdata~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",log="y",bty="L")
    points((Fns-baseline),col="green3",pty=16,cex=0.8)
    lines((Fns-baseline),col="green3")
    if(isTRUE(bootstrap)){
     #confidence interval
     lines(cycles,(CIFns[1,]-baseline),col="gold",lty=2);lines(cycles,(CIFns[2,]-baseline),col="gold",lty=2)
     #write some text in margin as legenda:
     mtext("original sample estimate",col="green3",cex=0.8,side=1,line=2,adj=-0.15)
     mtext("confidence limits",col="gold",cex=0.8,side=1,line=3,adj=-0.13)}   
    abline(h=ysdm) #plot threshold
    if(isTRUE(bootstrap)){
      #indicate Ct confidence interval
      lines(c(Ct,Ct),c(-1000,fluo[floor(Ct)]),lty=3,col="green3")
      lines(c(CICt[1],CICt[1]),c(-1000,CIFns[1,floor(CICt[1])]),lty=3,col="gold")
      lines(c(CICt[2],CICt[2]),c(-1000,CIFns[2,floor(CICt[2])]),lty=3,col="gold")}
    mtext("B",side=3,cex=1.5,adj=1)
    }
    
  ##Plot B: cycle-efficiency plots
  ######################################

  #open graphical device
    x11(width=14,height=7)
    par(mfcol=c(1,2))
  #Prepare for plotting    
    sdataplus <- sdata
    for(i in stada:1){sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],lp)))} 
    
  #depending on analysis print different plot:

  ##plots if plateau is reached:
  if(plotto==1){
    #optional plot 1
    plot(Lndata~sdata,ylim=c(-7,1),xlab="Fluorescence",ylab="ln(ln(Efficiency))",bty="L")
    if(isTRUE(bootstrap)){
      #confidence interval
      lines(sdata,log(log(CIEns[1,])),col="gold",lty=2);lines(sdata,log(log(CIEns[2,])),col="gold",lty=2)}    
    abline(v=sdata[stada],lty=3,col="blue");abline(v=sdata[stoda],lty=3,col="blue");abline(v=sdata[stp],lty=3,col="red")
    curve(Bob(x,bp),0,sp[3],add=TRUE,col="green3",lty=2)
    mtext("A",side=3,cex=1.5,adj=1)
    #abline(lp0);abline(lp1);abline(lp2)
        
    #optional plot 2
    plot(Endata~cycles,ylim=c(0.9,2.5),xlab="Cycle",ylab="Efficiency",bty="L")
    lines(cycles,Ens,lty=2,col="green3")
    if(isTRUE(bootstrap)){
     #confidence interval
     lines(cycles,CIEns[1,],col="gold",lty=2);lines(cycles,CIEns[2,],col="gold",lty=2)
     #write some text in margin as legenda:
     mtext("original sample estimate",col="green3",cex=0.8,side=1,line=2,adj=-0.15)
     mtext("confidence limits",col="gold",cex=0.8,side=1,line=3,adj=-0.13)}   
    mtext("B",side=3,cex=1.5,adj=1)
    }

    ##plots if plateau is NOT reached    
  if(plotto==2){
    #optional plot 1
    plot(Lndata~sdata,ylim=c(-7,1),xlab="Fluorescence",ylab="ln(ln(Efficiency))",bty="L")
    abline(v=sdata[stada],lty=3,col="blue");abline(v=sdata[stoda],lty=3,col="blue")
    curve(Bert(x,lp),0,sdata[stoda],add=TRUE,lty=2,col="green3")
    if(isTRUE(bootstrap)){
      #confidence interval
      lines(sdata[1:stoda],log(log(CIEns[1,1:stoda])),col="gold",lty=2);lines(sdata[1:stoda],log(log(CIEns[2,1:stoda])),col="gold",lty=2)}  
    mtext("A",side=3,cex=1.5,adj=1)
    
    #optional plot 2
    plot(Endata~cycles,ylim=c(0.9,2.5),xlab="Cycle",ylab="Efficiency",bty="L")
    lines(cycles,Ens,lty=2,col="green3")
    if(isTRUE(bootstrap)){
     #confidence interval
     lines(cycles[1:stoda],CIEns[1,1:stoda],col="gold",lty=2);lines(cycles[1:stoda],CIEns[2,1:stoda],col="gold",lty=2)
     #write some text in margin as legenda:
     mtext("original sample estimate",col="green3",cex=0.8,side=1,line=2,adj=-0.15)
     mtext("confidence limits",col="gold",cex=0.8,side=1,line=3,adj=-0.13)}   
    mtext("B",side=3,cex=1.5,adj=1)
    }    

}#end plotting section

########################################
##         Result Section             ##
########################################
    
    ############################################################################
    ## The output if all went well
    ############################################################################    
      #make sure baseline fits original data!: revert overbaseline substraction
      if(mark){sbase[1]<-sbase[1]-dummy} 
    #Actual output:
    if(isTRUE(bootstrap)){
      if(plotto==2){bp<-c(NA,NA,NA,NA,NA,NA)}#generate NA bp in case of truncation
      tabletop<-cbind(c(Emax,ai0,Ct,lp,bp),matrix(c(CIEmax,CIai0,CICt,as.vector(CIlip),as.vector(CIbip)),nrow=12,ncol=2,byrow=T))
      rownames(tabletop)<-c("Emax","ai0","Ct","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")
      colnames(tabletop)<-c("Estimate","lower limit","upper limit")
    return(tabletop)
    }else{
    if(output=="parameters"){
      if(plotto==1){return(list("baseline"=sbase,"5PLM"=sp,"Quadratic"=lp,"Bilinear"=bp))}
      if(plotto==2){return(list("baseline"=sbase,"5PLM"=sp,"Quadratic"=lp,"Bilinear"=rep(NA,times=6)))}}
    if(output=="estimates"){
      tzadaaam<-c(sdm,Emax,ai0)
      names(tzadaaam)<-c("Cq","Emax","a*i0")
      return(tzadaaam)}
    if(output=="all"){    
      if(plotto==1){tzadaaam<-c(sdm,Emax,ai0,sbase,sp,lp,bp)}
      if(plotto==2){tzadaaam<-c(sdm,Emax,ai0,sbase,sp,lp,rep(NA,times=6))}
      names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")    
      return(tzadaaam)}}
    
########################################
##          Error Section             ##
########################################

    #phase fit failed (no lp and/or bp)
    ########################
    }else{
    if(!isTRUE(silent)){
      message("Failed to obtain bilinear model")
      if(any(fail==301)){message("Linear model fit failed")}
      if(any(fail==302)){message("Parameter optimisation failed")}
      message("Please check data")
      message("----")}
    }#end of phase fitting control
    
    #5PLM fit failed
    ########################
    }else{   
    if(!isTRUE(silent)){
      if(all(fail!=204)){message("Failed to fit 5PLM")}
      if(any(fail==201)){message("Please check data")}
      if(any(fail==202)){message("Absence of amplification suspected")}
      if(any(fail==203)){message("Kinetic demarcation failed");message("Absence of amplification suspected")}
      if(any(fail==204)){message("5PLM residual error per cycle [resec value] exceeds tolerance");message("Absence of amplification suspected")}
      message("----")}
    }#end of 5PLM fit control 

    #Net Fluo increase too small!
    ########################
    }else{
    if(!isTRUE(silent)){
      message("Net fluorescence increase too low")
      message("Absence of amplification suspected")
      message("----")}
    }#end of fluo control if-else structure
    
    #Logical Syntax incorrect!
    ########################
    }else{
    if(!isTRUE(silent)){message("You're not making sense")
      if(!is.logical(plots)){message("'plots' is supposed to be of type 'logical'")}
      if(!is.logical(silent)){message("'silent' is supposed to be of type 'logical'");message("Printing this additional message just to annoy you...")}
      if(!is.logical(debugg)){message("'debugg' is supposed to be of type 'logical'")}
      message("----")}
    }#end of logical syntax control if-else structure
    
    #Output Syntax incorrect!
    ########################    
    }else{
    if(!isTRUE(silent)){message("You're not making sense")
      message("'output' should be either 'estimates', 'parameters', or 'all' ")
      message("----")}   
    }#end of baseline syntax control if-else structure

    #Baseline Syntax incorrect!
    ########################    
    }else{
    if(!isTRUE(silent)){message("You're not making sense")
      message("'baseline' should be either 'flat' or 'slanted' ")
      message("----")}   
    }#end of baseline syntax control if-else structure
    
    #NA output (default)!
    ########################  
    if(!isTRUE(silent)){message("We're sorry, there is no numerical output")
    if(isTRUE(plots)){message("We're sorry, but there are no plots either")}
    message("----")}
      sbase<-c(NA,NA);names(sbase)<-c("intercept","slope")
      sp<-c(NA,NA,NA,NA,NA,NA);names(sp)<-c("y0","fmax","xmid","b","g")
      lp<-c(NA,NA,NA);names(lp)<-c("int","slo1","slo2")     
      bp<-c(NA,NA,NA,NA,NA,NA);names(bp)<-c("a2","a1","a3","ipt","eta","chi")       
    if(isTRUE(bootstrap)){
    tabletop<- matrix(nrow=16,ncol=3,data=rep(NA,times=48))
    rownames(tabletop)<-c("Emax","ai0","Ct","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")
    colnames(tabletop)<-c("Estimate","lower limit","upper limit")
    return(tabletop)
    }else{
    if(output=="parameters"){
      return(list("baseline"=sbase,"5PLM"=sp,"Quadratic"=lp,"Bilinear"=bp))
    }else{
    if(output=="all"){    
      tzadaaam<-c(NA,NA,NA,sbase,sp,lp,bp)
      names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")    
      return(tzadaaam)
    }else{#default output
      tzadaaam<-c(NA,NA,NA)
      names(tzadaaam)<-c("Cq","Emax","a*i0")
      return(tzadaaam)}}}
      }#end of function
    
################################################################################    
###EXTRA! EXTRA!
    
##The following function is a shell
##its sole purpose is to prevent a function error from breaking multiple subsequent
##applications of 'analyse', eg. in the case of the application of the function to 
##a full array of PCR reactions: apply(data,2,analyse). Normally a function error would prevent
##all output, apply(data,2,semper) will ensure output of the results.

semper<-function(x,base.line="optimized",output="estimates",bootstrap=FALSE,plots=FALSE,silent=FALSE,kalman=TRUE,debugg=FALSE,n=1000,p=0.05,resec.tolerance=0.125){
        go<-try(analyse(x,base.line,output,bootstrap,plots,silent,kalman,debugg,n,p,resec.tolerance),silent=TRUE)
        if(inherits(go, "try-error")){
          sbase<-c(NA,NA);names(sbase)<-c("intercept","slope")
          sp<-c(NA,NA,NA,NA,NA,NA);names(sp)<-c("y0","fmax","xmid","b","g")
          lp<-c(NA,NA,NA);names(lp)<-c("int","slo1","slo2")     
          bp<-c(NA,NA,NA,NA,NA,NA);names(bp)<-c("a2","a1","a3","ipt","eta","chi")       
        if(isTRUE(bootstrap)){
          tabletop<- matrix(nrow=16,ncol=3,data=rep(NA,times=48))
          rownames(tabletop)<-c("Emax","ai0","Ct","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")
          colnames(tabletop)<-c("Estimate","lower limit","upper limit")
          return(tabletop)
        }else{
          if(output=="parameters"){
            return(list("baseline"=sbase,"5PLM"=sp,"Quadratic"=lp,"Bilinear"=bp))
          }else{
            if(output=="all"){    
              tzadaaam<-c(NA,NA,NA,sbase,sp,lp,bp)
              names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")    
              return(tzadaaam)
            }else{#default output
              tzadaaam<-c(NA,NA,NA)
              names(tzadaaam)<-c("Cq","Emax","a*i0")
              return(tzadaaam)}}}
        }else{
          return(go)}
        }

##The following function is a shell
##it has the same functionality as the above function but returns NaN in stead 
##of NA in case of failure of the FPK procedure. This makes it straight forward
## to distinguish between non-amplification events and actual function failure

semper2<-function(x,base.line="optimized",output="estimates",bootstrap=FALSE,plots=FALSE,silent=FALSE,kalman=TRUE,debugg=FALSE,n=1000,p=0.05,resec.tolerance=0.125){
        go<-try(analyse(x,base.line,output,bootstrap,plots,silent,kalman,debugg,n,p,resec.tolerance),silent=TRUE)
        if(inherits(go, "try-error")){
          sbase<-c(NaN,NaN);names(sbase)<-c("intercept","slope")
          sp<-c(NaN,NaN,NaN,NaN,NaN,NaN);names(sp)<-c("y0","fmax","xmid","b","g")
          lp<-c(NaN,NaN,NaN);names(lp)<-c("int","slo1","slo2")     
          bp<-c(NaN,NaN,NaN,NaN,NaN,NaN);names(bp)<-c("a2","a1","a3","ipt","eta","chi")       
        if(isTRUE(bootstrap)){
          tabletop<- matrix(nrow=16,ncol=3,data=rep(NaN,times=48))
          rownames(tabletop)<-c("Emax","ai0","Ct","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")
          colnames(tabletop)<-c("Estimate","lower limit","upper limit")
          return(tabletop)
        }else{
          if(output=="parameters"){
            return(list("baseline"=sbase,"5PLM"=sp,"Quadratic"=lp,"Bilinear"=bp))
          }else{
            if(output=="all"){    
              tzadaaam<-c(NaN,NaN,NaN,sbase,sp,lp,bp)
              names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","slo1","slo2","a2","a1","a3","ipt","eta","chi")    
              return(tzadaaam)
            }else{#default output
              tzadaaam<-c(NaN,NaN,NaN)
              names(tzadaaam)<-c("Cq","Emax","a*i0")
              return(tzadaaam)}}}
        }else{
          return(go)}
        }


##The following function provides use of both Brittish & American spelling of "analyse"
analyze<-function(x,baseline="slanted",output="estimates",bootstrap=FALSE,plots=FALSE,silent=FALSE,kalman=TRUE,debugg=FALSE,n=1000,p=0.05,resec.tolerance=0.125){
   return(analyse(x,baseline,output,bootstrap,plots,silent,kalman,debugg,n,p,resec.tolerance))} 
