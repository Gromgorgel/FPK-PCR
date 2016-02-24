## Disclaimer: This function has only been tested on a limited amount of data
##             there may be a few bugs left. Please contact the author with 
##             any remarks and/or suggestions.
##
## Original Publication: http://nar.oxfordjournals.org/content/40/2/e10.long
## Updates will appear at: github.com/Gromgorgel
##
## Most problems are caused by a lack of data resulting in a suboptimal fit
## try running more cycles and increasing primer concentration.
## If that doesn't help, feel free to contact me.

################################################################################
#The vertex-series of algorithms has the power to force a negative x-axis location
#of the vertex during the first phase fit. Vertex 2 only uses vertex control if
#standard fit failed.
#Version 11 has an updated sdataplus calculation that does not involve constant E
################################################################################

##LEGENDA
##plots = TRUE to make some nice output plots which may help locate problems
##baseline =  either "slanted" or "flat" depending on the model to be used
##output =  either "estimates" (for E, ai0 & Cq), "parameters" (for model parameters as a list), or "all" (both as a concatenated vector)
##silent= TRUE to prevent the algorithm from printing informative messages

analyse<-function(fluo,baseline="slanted",output="estimates",plots=FALSE,silent=FALSE){ 
    fluo<-as.numeric(fluo)                  #coerce numeric format
    endd<-length(fluo)                      #number of cycles ran
    cycles<-c(1:endd)                       #generate cycle numbers      
    plotto<-1                               #indicates which type of plots are generated (depends on available data, default=1)
    fail<-0                                 #keeps track of failures in order to return appropriate error message (default=0)
    vfail<-0
    mark<-FALSE                             #keeps track of over background subtraction
    require(minpack.lm)                     #load Levenberg-Marquardt nonlinear least-squares algorithm 
    require(grDevices)

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
if(baseline=="slanted"|baseline=="flat"){
#Output Syntax: check if user input is correct (output is either "estimates" , "parameters", or "all")
if(output=="estimates"|output=="parameters"|output=="all"){
#Logical Syntax: check if user input is correct (i.e. are the logicals indeed logicals)
if(is.logical(c(plots,silent))){

########################################
##            Data-Control            ##
########################################

#check if there is an obvious absence of amplification (i.e. net increase in 
#fluorescence is less than twofold, or mean fluorescence is negative)

#note that we use a rough baseline estimate & remove all NAs!
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
    if(!isTRUE(silent)){message("Missing values detected: some cycles have NA fluorescence");message("no action taken")}
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
  dummy<-2+abs(min(fluo,na.rm=T))   
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
  	a<-abs(min(fluo)-max(fluo))				  #calcute start parameter a
  	Y0<-mean(fluo[3:8])					        #calcute start parameter Y0
	  X0<-length(fluo[fluo<(Y0+a/2)])			#calcute start parameter X0
    ##special treatment for b
    bgrid <- expand.grid(X0= X0, Y0= Y0, a= a, b= seq(-20,-1,by=0.5))		 #define the grid to search best b in
    FPLM <- function(x){return(x[2]+x[3]/(1+(cycles/x[1])^x[4]))}		     #define FPLM to search b with
    billies <- apply(bgrid,1,FPLM)							                         #generate possible outcomes
    bestb <- which.min(colSums(apply(billies,2,function(p){p-fluo})^2))  #best outcome is that which minimizes the difference between real and predict
    b <- bgrid[bestb,4]
	  sp<-c(Y0= Y0,a= a,X0= X0,b= b)
    #We now try to optimize starting values by perfroming a 4PLM fit (this is crucial, bad initial parameter estimates will screw with the 5PLM fit)
    #more so since the nls.lm routine has the tendency to keep soldiering on and return crazy parameter values instead of simply giving up
    mleko <- deriv(~Y0+a/(1+(x/X0)^b),c("Y0","a","X0","b"),function(x,Y0,a,X0,b){})
	  modle<-try(nls(fluo~mleko(cycles,Y0,a,X0,b),start=sp,control=list(maxiter=100)),silent=TRUE)
	  #if the 4PLM fit fails, send message and use rough estimates
	  if(inherits(modle,"try-error")){
	       if(!isTRUE(silent)){message("Optimization of manual parameter values failed (4PLM fit)...")}
	       }else{
	       sp<-coef(modle)}		 
   #final message (separator)
   if(!isTRUE(silent)){message("----")}
  }
  sp<-c(sp[1],0.01,sp[2],sp[3],(0.335-0.0195*sp[4]),1)
  names(sp)<-c("y0","s","fmax","xmid","b","g")
##define model function
  fRichL<-function(x,p){return(p[3] + p[2]*x + (p[1] - p[3]) / (1 + (2^(1/p[6])-1) * exp(p[5]*(x-p[4])))^(p[6]))}
##define Residual Sum of Squares function for nls.lm
  sRRs<-function(m){resi<-fluo-fRichL(cycles,m)
                    return(resi)}
##fit 5PLM model to data
  sp<-try(nls.lm(par=sp,fn=sRRs),silent=TRUE)
  if(inherits(sp, "try-error")|any(is.nan(sp$par)==T)){
      sp<-c(NA,NA,NA,NA,NA,NA)
      fail<-201
      #Failure to fit a 5PLM may be a symptom of absence of amplification
      #if there is less than a 32 fold increase in fluorescence there is probably no amplification (less then 5 cycles of amplification)
      if((fluo[endd]/fluo[1])<32){fail<-202}                   
  }else{
  sp<-sp$par}
  names(sp)<-c("y0","s","fmax","xmid","b","g")


########################################
##          5PLM fit Check            ##
########################################

#We calculate the necessary values from the 5PLM fit

  #Second derivative maxima aka the Cq value(s)
    sdm<-as.numeric((sp[5]*sp[4]+log(-(1/2)*(-1-3*sp[6]+sqrt(1+6*sp[6]+5*sp[6]^2))/(sp[6]^2*(2^(1/sp[6])-1))))/sp[5])
    sdm2<-as.numeric((sp[5]*sp[4]+log((1/2)*(1+3*sp[6]+sqrt(1+6*sp[6]+5*sp[6]^2))/(sp[6]^2*(2^(1/sp[6])-1))))/sp[5])

  #calculate reference points of fluorescence accumulation
    fluospar<-function(p){return(floor(log(((-1/(p-1))^(1/sp[6])-1)/(2^(1/sp[6])-1))/sp[5]+sp[4]))}
    stada<-fluospar(0.05)+1           #cycle in which 5% of total fluorescence increase is reached
    stoda<-fluospar(0.85)             #cycle in which 85% of total fluorescence increase is reached
      #in case of data truncation the plateau normally follows directly after the truncation causing the 
      #95%FU mark to lie within the available cycle range, this may not always be the case
      #therefore it is useful to check if the 85% mark does not lie outside of the available cycle range
      if(stoda>endd){stoda<-endd}
    stp<-fluospar(0.95)               #cycle in which 95% of total fluorescence increase is reached
    startg<-4                         #start groundphase: skipping first three cycles to avoid high-noise data
    stopg <-stada-7                   #stop groundphase: skipping cycles that may be contaminated with amplicon fluorescence
    #return(c(startg,stopg,stada,stoda,stp))

#If any of these return NA that's a strong indication that the 
#data do not contain a valid PCR curve
if(any(is.na(c(stada,stoda,stp)))){fail<-203}

#if the fit of the 5PLM failed OR the values calcuated from the fit are invalid
#there is no use in contiuing the procedure. Therefore we only proceed if fail equals zero
if(fail==0){

########################################
##        Baseline Analysis           ##
########################################

##check for very "early" reactions (Cq<10, too few cycles for baseline fitting)
  if(baseline=="flat"|(stopg-startg)<1){
    if(!isTRUE(silent)&(stopg-startg)<1){message("Groundphase too short for linear model!");message("Using flat baseline instead...");message("----")}    
    sbase<-c(mean(fluo[startg:stopg]),0)
    }else{
##fit baseline (linear model: least squares, normal error).
    sbase<-lm(fluo[startg:stopg]~cycles[startg:stopg])$coefficients}
    names(sbase)<-c("intercept","slope")
    baseline<-sbase[2]*cycles+sbase[1]
    sdata<-fluo-baseline 
    #in Phase 2 the double log often causes NaNs which we cannot use in therefore we resort to
    #using abs() as a means of data imputation to counteract loss of points due to NaN from log(-)   
    Lndata<-c(NA,log(abs(log(sdata[-1]/sdata[-endd]))))
################################################################################     
                ########################################
                ##         End of Pre-Analysis        ##
                ##    Start of Calculation Section    ##
                ########################################
################################################################################     

##General Function Repository:

  #Define a function that calculates the vertex and thus allows checking its x-axis location
     vertex<-function(lp){output<-c(-lp[2]/(2*lp[3]),(4*lp[3]*lp[1]-lp[2]^2)/(4*lp[3]))
                         names(output)<-c("h","k")
                         return(output)} 
  #Define a function that can swap model paramters from vertical to horizontal and back
  swapper<-function(qp){
              if(length(qp)==2){ #swap from horizonatal to vertical:
                aa<--1/qp[2]^2; bb<--2*aa*qp[1]; cc<-(1/4)*(bb^2)/aa; return(c(cc,bb,aa))}
              if(length(qp)==3){ #swap from vertical to horizontal:
                aa<--qp[2]/(2*qp[3]);bb<-sqrt(-1/qp[3]);return(c(aa,bb))}}
  #Define Mono-linear model 
    Bert<-function(n,l){return( l[1]+n*l[2]+l[3]*n^2 )}   
  #Define the bilinear model
     Bob<-function(n,b){return(b[6]+b[5]*log(exp((b[1]*(n-b[4])^2+b[2]*(n-b[4]))/b[5])+exp(b[3]*(n-b[4])/b[5])))}                          


##Irrespective of the number of phases in the data, the fit of the first phase is always necessary and identical in process:

        #select data cycles & respective fluorescence values
        Fnsub<-sdata[stada:stoda]
        Lnsub<-Lndata[stada:stoda]
        
## some notes on what is about to follow
## firstly, a simple linear (quadratic) model is fit to the data (lp0) BUT two things are not entirely right:
## (1) when fitting the linear model to the first E decline phase it makes more sense to
##     minimize the fluorescence residuals (horizontal distance) since cycle efficiencies are 
##     themselves calculated from the fluorescence measurements.           
## (2) There is no guarantee that the vertex has a negative x-axis location. if the vertex has
##     a positive location AND the curve is convex than the efficiency rises before it declines. That's just silly.
##     Therefore we: (1) try minimizing the vertical distances, but when that fails we try minimizing the 
##     horizontal distances instead. (2) IF the vertex of either if (horiz or vertic) is x-positive (and the curve is convex) we fit 
##     a modified form of the quadratic equations which ensures the vertex has a negative x-axis location 

##Start of the fitting story: we try to minimize the horizontal distances
      lpM<-try(lm(Fnsub~I(abs(Lnsub)^(1/2))),silent=TRUE)
      #If the horizontal distances didn't work we try the vertical distances
      if(inherits(lpM, "try-error")){
         fail<-append(fail,303)
         #use vertical distances
         lp0<-try(coef(lm(Lnsub~Fnsub+I(Fnsub^2))),silent=T)
         #If the lp0 didn't work the story ends here
         if(inherits(lp0,"try-error")){      
          fail<-append(fail,301)       
         #If the lp0 fit worked we continue by checking the vertex
         }else{         
          #if the vertex already has a negative location there is no need to recaculate, that would yield the same result anyway
          #so a quick check to see if the vertex indeed has a negative x-axis location and if the curve is convex (i.e. quadratic term is negative)
          #if not we re-fit using the adjusted formulae:
          if((vertex(lp0)[1]>0)&(lp0[3]<0)){                #IF vertex has a positive x-axis location, we force a negative one:
                              vertfor<-function(n,p){return(p[1]*(n+p[2]^2)^2+p[3])}  #this is the forced negative vertex formulae
                              vRSS<-function(x){return(Lnsub-vertfor(Fnsub,x))}       #this function returns the residuals
                              vstar<-c(lp0[3],sqrt(abs(lp0[2]/(2*lp0[3]))),(4*lp0[3]*lp0[1]-lp0[2]^2)/(4*lp0[3])) #these are the starting parameter estimates (based on lp0 fit)
                              #actual fit: if it fails we stick to the lp0 we found earlier & remember the failure to give an accurate message
                              verti<-try(nls.lm(fn=vRSS,par=vstar),silent=TRUE)               
                              if(inherits(verti, "try-error")|verti$info==9){ #error 9= maxiter reached
                                  fail<-append(fail,305)
                                }else{
                                  verti<-coef(verti)
                                  lp0<-c((verti[3]*4*verti[1]+(-(2*verti[1]*verti[2]^2))^2)/(4*verti[1]),-(2*verti[1]*verti[2]^2),verti[1])}
                              }#and of lp0 vertex control
          }#end of lp0 vertex check
          lp<-lp0          
      #If the horizontal distances worked we continue by converting the parameters and checking the vertex          
      }else{                
        #conversion of parameters vertical distance space:
        lp<-swapper(coef(lpM))
        
        #if the vertex already has a negative location there is no need to recaculate, that would yield the same result anyway
        #so a quick check to see if the vertex indeed has a negative x-axis location and negative quadratic term (convex)
        #if not we re-fit using the adjusted formulae
        if((vertex(lp)[1]>0)&(lp[3]<0)){
                Vertfor<-function(n,p){return(sqrt((n-p[3])/p[1])+p[2]^2)}  #this is the forced negative vertex formulae
                VRSS<-function(x){return(Fnsub-Vertfor(Lnsub,x))}           #this function returns the residuals
                Vstar<-c(lp[3],sqrt(abs(lp[2]/(2*lp[3]))),(4*lp[3]*lp[1]-lp[2]^2)/(4*lp[3])) #these are the starting parameter estimates (based on lp0 fit)
                Verti<-try(nls.lm(fn=VRSS,par=Vstar),silent=TRUE)
                if(inherits(Verti, "try-error")|Verti$info==9){
                      fail<-append(fail,305)
                }else{
                      Verti<-coef(Verti)
                      lp<-c((Verti[3]*4*Verti[1]+(-(2*Verti[1]*Verti[2]^2))^2)/(4*Verti[1]),-(2*Verti[1]*Verti[2]^2),Verti[1])
                }#and of lpM vertex control
        }#end of lpM vertex check
     }#END of the fitting story we either have Fail or lp

################################################################################
# we should only continue of we have parameters (lp exists, i.e. the first phase fit succeeded)!
if(exists("lp")){
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
    lp2<-coef(lm(Fnsub~Lnsub))
      #the probability exists that the linear model wil yield a positive slope
      #which causes a contradiction in the bilinear model
      #in those cases we force a negative slope through a custom fitting approach
      if(lp2[2]>=0){
      B.negative<-function(zz){return(sum((Fnsub-(zz[1]-zz[2]^2*Lnsub))^2,na.rm=TRUE))}
      lp2<-optim(par=c(Fnsub[1],0),fn=B.negative)$par
      lp2[2]<--lp2[2]^2}
    #convert the parameters for so the plot remains E vs Fn
    lp2<-c((-lp2[1]/lp2[2]),(1/lp2[2]))  
        
## Second: construct the bilinear model using the linear fits from both phases             
################################################################################     

  #estimate where both phases cross 
    crx <- (1/2)*(lp2[2]-lp[2]+sqrt(lp2[2]^2-2*lp2[2]*lp[2]+lp[2]^2+4*lp[3]*lp2[1]-4*lp[3]*lp[1]))/lp[3] 
    
  ##we now gather all intial parameters for the bilinear model
    slo1a<-lp[3]
    slo1b<-lp[2]+2*crx*lp[3]    #see remark:
        #for simple linear components of the bilinear model, the slope is not affected by the (x-crx) term (displacement of the function along the x-axis).
        #Only the intercept changes when the function is displaced by 'crx', and the intercept is not used in the bilinear model.
        #For quadratic equations the second slope (b) and intercept (int) change by displacement along the x-axis: int + b*(x+crx) + a*(x+crx)^2 <=> (int+b*crx+a*crx^2) + x * (b+2*crx*a) + a*x^2
        #thus displacement of the curve by 'crx' changes the 'b' slope, hence we adjust its value.
    slo2<-lp2[2]  
    ipt<- crx
    eta<--4
    chi<-lp[1]-eta*log(exp((slo1a*ipt^2-slo1b*ipt)/eta)+exp(-slo2*ipt/eta))
    bp<-c(slo1a,slo1b,slo2,ipt,eta,chi)

##eta and crx can do with a little optimizing
    #redifine action zone
    Fnsub<-sdata[stoda:endd]
    Lnsub<-Lndata[stoda:endd]      

    bobslee<-function(d){
      return(
      sum((Lnsub-
         (lp[1]-d[2]*log(exp((bp[1]*d[1]^2-bp[2]*d[1])/d[2])+exp(-bp[3]*d[1]/d[2]))
         +d[2]*log(exp((bp[1]*(Fnsub-d[1])^2+bp[2]*(Fnsub-d[1]))/d[2])+exp(bp[3]*(Fnsub-d[1])/d[2])))      
          )^2,na.rm=TRUE)
          )}
    crash<-try(optim(fn=bobslee,par=bp[4:5],method="L-BFGS-B",lower=c(max(sdata)*0.8,-8),upper=c(max(fluo),-3)),silent=TRUE)
    #this optimalisation may crash, in that case we keep the original values and take notice of the error
    if(inherits(crash, "try-error")){
    fail<-append(fail,307)
    }else{
    bp[4:5]<-crash$par
    #update slo1b
    bp[2]<-lp[2]+2*bp[4]*lp[3]
    #update Chi
    bp[6]<-lp[1]-bp[5]*log(exp((bp[1]*bp[4]^2-bp[2]*bp[4])/bp[5])+exp(-bp[3]*bp[4]/bp[5]))}
    names(bp)<-c("a2","a1","a3","ipt","eta","chi")   
   
##we end this section by defining a function for the bilinear model and calculating its intial Efficiency value   
    EmaxB0<-exp(exp(lp[1]))

## Third: define the fucntion to rebuild the fluorescence readings using the calculated cycle efficiencies  (TWO phase version)
################################################################################     

##Here comes the Residual sum of squares function
##it takes the 'fluorescence due to initial target copies' as input value since that is the value we want to optimize
    mini_one<-function(m){return(
        #we'll calculate the sum of squared residuals for the holistic model, so we'll start with the sum command 
				sum(weight*                                                                    
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
          cumprod(exp(exp(Bob(sdataplus,bp))))                                                            
        #end of cycle Fluorescence calculation
        )                                                        
				#now we just have to sqaure the residuals and sum them up
        )^2)
    ##that's all folks
    )}
    

################################################################################
#In case no second phase could be detected:
}else{ #part of the data truncation structure 
################################################################################
  #Inform user that data are truncated and no plateau was detected
    if(!isTRUE(silent)){message("Failed to detect second phase");message("Analysis limited to first phase of efficiency decline...");message("----")}
  
  #write a small function that describes the first phase
    EmaxB0<-exp(exp(lp[1]))
  #we also change the plotting output, so something useful can be produced
    plotto<-2

## define the fucntion to rebuild the fluorescence readings using the calculated cycle efficiencies  (ONE phase version)
################################################################################     

  ##the incomplete data case requires a dedicated Residual sum of squares function
  ##it takes the 'fluorescence due to initial target copies' as input value since that is the value we want to optimize.
  ##Also: it does not use the data beyond the 85% fluorescence mark since the since validity of the efficiency model
  ##beyond that point is not guaranteed
    mini_one<-function(m){return(
        #we'll calculate the sum of squared residuals for the holistic model, so we'll start with the sum command
				sum(weight[1:stoda]*
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
          cumprod(exp(exp(Bert(sdataplus[1:stoda],lp))))
        #end of cycle Fluorescence calculation
        )
				#now we just have to sqaure the residuals and sum them up
        )^2)
    ##that's all folks
    )}
    
################################################################################    
}#END of the data truncation structure (from here on, all commands are valid irrespective of the number of phases in the data)
################################################################################

########################################
##        Estimation of a*i0          ##
########################################

##In the baseline subtracted data we replace the groundphase by an estimate of the baseline subtracted fluo based on the Emax we found
##this prevents the noise in the baseline from affecting the ai0 estimate
  #first copy baseline subtracted data
  sdataplus <- sdata
  #now replace groundphase cycles with kinetic estimates
  for(i in stada:2){
    sdataplus[i-1]<-sdataplus[i]/exp(exp(Bert(sdataplus[i],lp)))} #Bert will do for both mono & bi linear since only the very beginning of the first phase is involved   
#define weights for ai0 determination process
    weight<-1.0001-sdataplus/sdataplus[endd]                                    
##The upper and lower limit for the initial target copies * amplicon Fluorescence are 1 and 1e-20 respectivly
    opti<-optimize(mini_one,interval=c(1e-10,1e+10))
    #we extract residual sum of squares, this allows calculation of the AIC
    #RSS<-opti$objective
    #aic<-endd+log(2*pi)+endd*log(RSS/endd)+2*(8+1)
    #we extract the actual result (value of ai0 that minimizes the sum of squares)
    result_one<-opti$minimum	   
    #we multiply with the scaling factor for the final result
    ai0<-result_one*1e-10

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
    ##Calculate the datapoints
     En<-(exp(exp(Bob(sdataplus,bp))))
     Fns<-baseline +                                                                  
          result_one *  1e-10 *
          cumprod(exp(exp(Bob(sdataplus,bp))))                                                        
    smother<-splinefun(cycles,sdata,method="fmm")
    threshold<-smother(sdm)

    #optional plot 1: fitted model on normal scale
    plot(fluo~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",bty="L")
    points(Fns,col="green3",pty=16,cex=0.8)
    lines(Fns,col="green3")
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
    abline(h=threshold)
    mtext("B",side=3,cex=1.5,adj=1)
    }
    
  ##plots if plateau is NOT reached    
  if(plotto==2){
    ##Calculate the datapoints
     En<-(exp(exp(Bert(sdataplus,lp))))
     Fns<-baseline +                                                                  
          result_one *  1e-10 *
          cumprod(exp(exp(Bert(sdataplus,lp)))) 
    smother<-splinefun(cycles,sdata,method="fmm")
    threshold<-smother(sdm)

    #optional plot 1: fitted model on normal scale
    plot(fluo~cycles,main="data and fitted model",xlab="cycles",ylab="fluorescence",bty="L")
    points(Fns,col="green3",pty=16,cex=0.8)
    lines(Fns,col="green3")
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
    abline(h=threshold)
    mtext("B",side=3,cex=1.5,adj=1)
    }
    
  ##Plot B: cycle-efficiency plots
  ######################################

  #open graphical device
    x11(width=14,height=7)
    par(mfcol=c(1,2))
  #Prepare for plotting    
    trat<-(sdata[2:endd])
    trut<-sdata[2:endd]/sdata[1:(endd-1)] 
  #depending on analysis print different plot:

  ##plots if plateau is reached:
  if(plotto==1){
    #optional plot 1
    plot(log(log(trut))~trat,ylim=c(-7,1),xlab="Fluorescence",ylab="ln(ln(Efficiency))",bty="L")
    abline(v=sdata[stada],lty=3,col="blue")
    abline(v=sdata[stoda],lty=3,col="blue")
    abline(v=sdata[stp],lty=3,col="red")
    curve(Bob(x,bp),0,sp[3],add=TRUE,col="green3",lty=2)
    mtext("A",side=3,cex=1.5,adj=1)
    
    #optional plot 2
    Ees<-exp(exp(Bob(sdataplus,bp)))
    plot(trut~cycles[-1],ylim=c(0.9,2.5),xlab="Cycle",ylab="Efficiency",bty="L")
    lines(cycles,Ees,lty=2,col="green3")
    mtext("B",side=3,cex=1.5,adj=1)
    }

    ##plots if plateau is NOT reached    
  if(plotto==2){
    #optional plot 1
    plot(log(log(trut))~trat,ylim=c(-7,1),xlab="Fluorescence",ylab="ln(ln(Efficiency))",bty="L")
    abline(v=sdata[stada],lty=3,col="blue")
    abline(v=sdata[stoda],lty=3,col="blue")
    curve(Bert(x,lp),0,sdata[stoda],add=TRUE,lty=2,col="green3")
    mtext("A",side=3,cex=1.5,adj=1)
    
    #optional plot 2
    Ees<-exp(exp(Bert(sdataplus,lp)))    
    plot(trut~cycles[-1],ylim=c(0.9,2.5),xlab="Cycle",ylab="Efficiency",bty="L")
    lines(cycles,Ees,lty=2,col="green3")
    mtext("B",side=3,cex=1.5,adj=1)
    }    

}#end plotting section

########################################
##         Result Section             ##
########################################
    
  ##print general warning messages if appropriate
    if(!isTRUE(silent)){
    if(any(fail==305)){message("Vertex master control failed");message("releasing constraints...")}
    if(any(fail==303)){message("Use of horizontal distances failed on 1st E decline phase");message("Minimizing vertical distances instead...")}
    if(any(fail==307)){message("Optimalistion of 'ipt' and 'chi' failed");message("Using non optimized starting values instead...")}
    message("----")}
    #this message will appear even if Silent=T (when appropriate)
    if(vertex(lp)[1]>0){message("WARNING: Vertex has positive x-axis location!")}

    ############################################################################
    ## The output if all went well
    ############################################################################    
      #make sure baseline fits original data!: revert overbaseline substraction
      if(mark){sbase[1]<-sbase[1]-dummy} 
    #Actual output:
    if(output=="parameters"){
      if(plotto==1){return(list("baseline"=sbase,"5PLM"=sp,"Phase1"=lp,"Bilinear"=bp))}
      if(plotto==2){return(list("baseline"=sbase,"5PLM"=sp,"Phase1"=lp,"Bilinear"=rep(NA,times=6)))}}
    if(output=="estimates"){
      tzadaaam<-c(sdm,EmaxB0,ai0)
      names(tzadaaam)<-c("Cq","Emax","a*i0")
      return(tzadaaam)}
    if(output=="all"){    
      if(plotto==1){tzadaaam<-c(sdm,EmaxB0,ai0,sbase,sp,lp,bp)}
      if(plotto==2){tzadaaam<-c(sdm,EmaxB0,ai0,sbase,sp,lp,rep(NA,times=6))}
      names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","sloA","sloB","a2","a1","a3","Fc","eta","Chi")    
      return(tzadaaam)}
    
########################################
##          Error Section             ##
########################################

    #phase fit failed (no lp)
    ########################
    }else{
    if(!isTRUE(silent)){
      message("Failed to fit bilinear model")
      message("Please check data")
      message("----")}
    }#end of phase fitting control
    
    #5PLM fit failed
    ########################
    }else{   
    if(!isTRUE(silent)){
      message("Failed to fit 5PLM")
      if(fail==201){message("Please check data")}
      if(fail==202){message("Absence of amplification suspected")}
      if(fail==203){message("Kinetic demarcation failed");message("Absence of amplification suspected")}
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
      if(!is.logical(silent)){message("'silent' is supposed to be of type 'logical'");message("Printing this extra message just to annoy you...")}
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
      lp<-c(NA,NA,NA);names(lp)<-c("int","sloA","sloB")     
      sp<-c(NA,NA,NA,NA,NA,NA);names(sp)<-c("y0","fmax","xmid","b","g")
      bp<-c(NA,NA,NA,NA,NA,NA);names(bp)<-c("a2","a1","a3","Fc","eta","chi")       
    if(output=="parameters"){
      return(list("baseline"=sbase,"5PLM"=sp,"Phase1"=lp,"Bilinear"=bp))
    }else{
    if(output=="all"){    
      tzadaaam<-c(NA,NA,NA,sbase,sp,lp,bp)
      names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","sloA","sloB","a2","a1","a3","Fc","eta","Chi")    
      return(tzadaaam)
    }else{#default output
      tzadaaam<-c(NA,NA,NA)
      names(tzadaaam)<-c("Cq","Emax","a*i0")
      return(tzadaaam)}}
      }#end of function
    
################################################################################    
###EXTRA! EXTRA!
    
##The following function is a shell
##its sole purpose is to prevent a function error from breaking multiple subsequent
##applications of 'analyse', eg. in the case of the application of the function to 
##a full array of PCR reactions: apply(data,2,analyse). Normally a function error would prevent
##all output, apply(data,2,semper) will ensure output of the results.

semper<-function(x,baseline="slanted",output="estimates",plots=FALSE,silent=FALSE){
        go<-try(analyse(x,baseline,output,plots,silent),silent=TRUE)
        if(inherits(go, "try-error")){
          sbase<-c(NA,NA);names(sbase)<-c("intercept","slope")
          lp<-c(NA,NA,NA);names(lp)<-c("int","sloA","sloB")     
          sp<-c(NA,NA,NA,NA,NA,NA);names(sp)<-c("y0","fmax","xmid","b","g")
          bp<-c(NA,NA,NA,NA,NA,NA);names(bp)<-c("a2","a1","a3","Fc","eta","chi")       
          if(output=="parameters"){
            return(list("baseline"=sbase,"5PLM"=sp,"Phase1"=lp,"Bilinear"=bp))
          }else{
            if(output=="all"){    
              tzadaaam<-c(NA,NA,NA,sbase,sp,lp,bp)
              names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","sloA","sloB","a2","a1","a3","Fc","eta","Chi")    
              return(tzadaaam)
            }else{#default output
              tzadaaam<-c(NA,NA,NA)
              names(tzadaaam)<-c("Cq","Emax","a*i0")
              return(tzadaaam)}}
        }else{
        return(go)}
        }

##The following function is a shell
##it has the same functionality as the above function but returns NaN in stead 
##of NA in case of failure of the FPK procedure. This makes it straight forward
## to distinguish between non-amplification events and actual function failure

semper2<-function(x,baseline="slanted",output="estimates",plots=FALSE,silent=FALSE){
        go<-try(analyse(x,baseline,output,plots,silent),silent=TRUE)
        if(inherits(go, "try-error")){
          sbase<-c(NaN,NaN);NaNmes(sbase)<-c("intercept","slope")
          lp<-c(NaN,NaN,NaN);NaNmes(lp)<-c("int","sloA","sloB")     
          sp<-c(NaN,NaN,NaN,NaN,NaN,NaN);NaNmes(sp)<-c("y0","fmax","xmid","b","g")
          bp<-c(NaN,NaN,NaN,NaN,NaN,NaN);NaNmes(bp)<-c("a2","a1","a3","Fc","eta","chi")       
          if(output=="parameters"){
            return(list("baseline"=sbase,"5PLM"=sp,"Phase1"=lp,"Bilinear"=bp))
          }else{
            if(output=="all"){    
              tzadaaam<-c(NaN,NaN,NaN,sbase,sp,lp,bp)
              names(tzadaaam)<-c("Cq","Emax","a*i0","intercept","slope","y0","s","Fmax","xmid","b","g","int","sloA","sloB","a2","a1","a3","Fc","eta","Chi")    
              return(tzadaaam)
            }else{#default output
              tzadaaam<-c(NaN,NaN,NaN)
              names(tzadaaam)<-c("Cq","Emax","a*i0")
              return(tzadaaam)}}
        }else{
        return(go)}
        }

##The following function provides use of both Brittish & American spelling of "analyse"
analyze<-function(x,baseline="slanted",output="estimates",plots=FALSE,silent=FALSE){return(analyse(x,baseline,output,plots,silent))}
