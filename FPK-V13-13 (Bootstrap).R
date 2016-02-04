#######Bootstrap procedure based on FFF-V08

# Fn            Fluo values of the data to be boostrapped
# Ln            log(log(En)) values of the data to be bootstrapped
# RSS           Residual sum of squares of the original sample linear model, needed to set the frequency cutoff
# KAL           T or F (Kalman Filter ON or OFF)

# Bpar1 & Jpar1 quadratic model parameters of the first phase boostraps & Jackknives, needed to calculate the crx-point(if block=2)

# ose           Original sample esitmates c(lp,lp2,bp), required to calculate the full set of Jackknives for the trilinar model(if block=2)
# block         either 1 or 2, sets bootstrapping to first or second phase of decline

# n             number of bootstraps performed
# silent        set to TRUE to turn of messages


phaseboot<-function(Fn,Ln,dmc=NA,RSS,KAL=FALSE,Evar=NA,Bpar1=NA,Jpar1=NA,ose=NA,lice=NA,block=1,n=1000,silent=FALSE){
 #some stuff we need:
 Fs  <-1/60*1000                      #Sampling Frequency: 1 measurement every 60 seconds, expressed in mHz
################################################################################     
                ########################################
                ##       Start of the Bootstrap's     ##
                ##     General Function Repository    ##
                ########################################
################################################################################     
 
  #Restricted freedom model's derivative
    mleko <- deriv(~-(a/b)^2-(a/b^2)*x-(2/b^2)*x^2,c("a","b"),function(x,a,b){})
  #The quadratic model 
    Bert<-function(n,l){return( l[1]+n*l[2]+l[3]*n^2 )}   
  #The bilinear model
    Bob<-function(n,b){return(b[6]+b[5]*log(exp((b[1]*(n-b[4])^2+b[2]*(n-b[4]))/b[5])+exp(b[3]*(n-b[4])/b[5])))}   

  #The Fourier transform PCR-lowpass filter
    lowpass<-function(y,cutoff,RSoS=F){
      #ensure circular continuity
      temp<-c(y[floor(length(y)/2):1],y,y[length(y):floor(length(y)/2+1)])
      tren<-coef(lm(temp[c(1,length(temp))]~c(1,length(temp))))  #calculate trendline
      temp<-temp-(tren[1]+c(1:length(temp))*tren[2])             #subtract trendline
      #Padding to get length 2^n (if necessary)
      if(!any((length(temp)/2^c(1:7))==1)){
        L<-(2^which((length(temp)/2^c(1:7))<1)[1])-length(temp) #number of zeroes needed
        temp<-append(temp,rep(0,times=L))}
      #Apply smoothing
                 lefou<-fft(temp) 
                #Frequency bins
                  freqs<-c((length(temp)/2):1)*Fs/(length(temp))     
                #find cutoff
                  cut.point<-tail(which(freqs>cutoff),1)   
                #check if anything has to ber removed 
                  not<-c((length(lefou)/2-cut.point+1):(length(lefou)/2+cut.point))
                #delete info beyond cutoff
                  lefou[not]<-0+0i
                #inverse transform
                  temp<-Re(fft(lefou,inverse=T)/length(lefou))
      #select smoothed data (fitted values)
      fits<-(temp+(tren[1]+c(1:length(temp))*tren[2]))[floor(length(y)/2+1):floor(3*length(y)/2)]
      sres<-sum((y-fits)^2)
      #output
      if(isTRUE(RSoS)){return(sres)}else{return(fits)}}
    
  #A function that generates the Hat-matrix (or Smoothing Matrix) for our smoother
    mad.hatter<-function(fop){
            #generate matrix of appropriate dimensions
            ttt<-length(fop)
            In <- diag(ttt) ## identity matrix
            Snw <- matrix(0, nrow = ttt, ncol = ttt) ##empty hat matrix
            #hat matrix (S) of a linear smoothing algorithm can be obtained by smoothing unit vectors
            #the result of smoothing the ith unit vector is the ith column of S
             for(j in 1:ttt){
              y <- In[,j]
              Snw[,j] <- lowpass(y, cutfreq)}
             #output
             return(Snw)}
             
 #A function to fit the first phasre to a Ln data vetcor (dator)
    ph1fit<-function(dator){
               Fsub<-Fn[dmc[3]:dmc[4]]
               Lsub<-Ln[1:(dmc[3]-1)]
               tlpm<-try(lm(dator~Fsub+I(Fsub^2)),silent=T)#fit quadratic model
               if(inherits(tlpm,"try-error")){#if FAIL we stop here
                pp<-c(NA,NA,NA)
               }else{# if SUCCESS we continue  by fitting the restricted modle
                #extract & transform the model parameters in to restricted parameter estimates
                pp<-coef(tlpm)
                lap<-c(sqrt(abs(pp[1]))*sqrt(1/abs(pp[3])),sqrt(1/abs(pp[3]))) ;  names(lap)<-c("a","b")
               #Restricted fit
                todle<-try(nls(dator~mleko(Fsub,a,b),start=lap,control=list(maxiter=100)),silent=TRUE)
                if(!inherits(todle,"try-error")){#if NO FAIL we continue (else the unrestriced pp estimate will be returned)
                  lap<-coef(todle)
                  pp<-c(-(lap[1]/lap[2])^2,-(lap[1]/lap[2]^2),-(2/lap[2]^2))
                  if(isTRUE(KAL)){#apply KALMfAN filter if necessary
                  ##We calculate fluorescence increases between each cycle
                    dFn<-c(NA,abs(Fn[1:(dmc[3]-1)]-Fn[2:dmc[3]]))       
                  ##Calculate the initial states for the filter to start from
                    #current state
                    currstate<-matrix(c(Bert(Fn[dmc[3]],pp),pp[2]+2*pp[3]*Fn[dmc[3]]),nrow=2,ncol=1)
                    states<-currstate
                    #we can calculate sd of the acceleration a using the residuals from the regression  (matrix form: sigma^2 =  MSE * (X'X)^-1, see page 207 of book)
                      xix<-matrix(c(rep(1,times=length(Fn[])),Fn,Fn^2),nrow=length(Fn),ncol=3,byrow=F)
                      MSE<-sum((straps[dmc[3]:dmc[4],a]-Bert(Fsub,pp))^2)/(length(Fsub)-3)   #we estimted three parameters so we loose 3 df
                      sig<- 2 * (MSE*diag(solve(t(xix)%*%xix)))[3]^(1/2)     # * 2 since our parameter c equals 2 a (acceleration) in the physical model
                    #variance-covariance matrix of the estimates
                    estcov<-matrix(c((sig*dFn[dmc[3]]^2)^2,sig^2*2*dFn[dmc[3]]^3,sig^2*2*dFn[dmc[3]]^3,(sig*2*dFn[dmc[3]])^2),nrow=2,ncol=2,byrow=T)
                  ##Actual Filter
                    for(i in dmc[3]:2){
                    ##coefficient matrices
                      AA<-matrix(c(1,-dFn[i],0,1),nrow=2,ncol=2,byrow=T)
                      BB<-matrix(c(dFn[i]^2,-2*dFn[i]),nrow=2,ncol=1)
                      CC<-matrix(c(1,0),nrow=1,ncol=2)
                    ##Covariance matrices
                      EE<-matrix(c((sig*dFn[i]^2)^2,sig^2*2*dFn[i]^3,sig^2*2*dFn[i]^3,(sig*2*dFn[i])^2),nrow=2,ncol=2,byrow=T)
                      EZ<-Evar[i]
                    ##The run
                      statepred<-AA%*%currstate+BB%*%pp[3]
                      measupred<-CC%*%statepred
                      estcov<-AA%*%estcov%*%t(AA)+EE
                      K<-estcov%*%t(CC)%*%solve((CC%*%estcov%*%t(CC)+EZ))
                      statepred<-statepred+K%*%(Lsub[i-1]-measupred)
                      states<-cbind(states,statepred)
                      estcov<-(diag(2)-K%*%CC)%*%estcov
                    ##refit model
                      Fnt<-Fn[(i-1):dmc[4]]
                      Lnt<-c(rev(states[1,]),Ln[(dmc[3]+1):dmc[4]])
                  	  todle<-nls(Lnt~mleko(Fnt,a,b),start=lap,control=list(maxiter=100))
                      lpK<-coef(todle)
                      lpK<-c(-(lpK[1]/lpK[2])^2,-(lpK[1]/lpK[2]^2),-(2/lpK[2]^2))
                    ##recalculate sig
                      xix<-matrix(c(rep(1,times=length(Fnt)),Fnt,Fnt^2),nrow=length(Fnt),ncol=3,byrow=F)
                      MSE<-sum((Lnt-Bert(Fnt,lpK))^2)/(length(Fnt)-3)
                      sig<- 2 * (MSE*diag(solve(t(xix)%*%xix)))[3]^(1/2)
                    ##advance current state using updated model
                      currstate<-matrix(c(Bert(Fn[i-1],lpK),lpK[2]+2*lpK[3]*Fn[i-1]),nrow=2,ncol=1)
                    }#END of the filter for-loop
                    #we check if the last Kalman model is valid and, if so, update the model paramters
                    if(!any(is.na(lpK))){pp<-lpK}
              }}}#close all IF/ELSE structeres
              return(pp)}
                                  
 #A functions to assemble the bilinear model
    bpmaker<-function(mp,mp2){  
             crx <- (1/2)*(mp2[2]-mp[2]+sqrt(mp2[2]^2-2*mp2[2]*mp[2]+mp[2]^2+4*mp[3]*mp2[1]-4*mp[3]*mp[1]))/mp[3] 
             ##finally, we gather all intial parameters for the bilinear model
             slo1a<-mp[3]
             slo1b<-mp[2]+2*crx*mp[3]    #see remark:
             slo2<-mp2[2]  
             eta<-rnorm(1,mean=ose[10],sd=0.1)  #since we are unable to optimize eta we just use a random nr around the ose to at least give some confidence interval
             chi<-mp[1]-eta*log(exp((slo1a*crx^2-slo1b*crx)/eta)+exp(-slo2*crx/eta))             
             bbp<-c(slo1a,slo1b,slo2,crx,eta,chi)             
             names(bbp)<-c("a2","a1","a3","ipt","eta","chi")                
             return(bbp)}

################################################################################
                ########################################
                ##      End of Function repository    ##
                ##      Start of Bootstrap section    ##
                ########################################
################################################################################

#Prepping: make all fitted values & residuals
#############################################
if(block==1){ #expects input of all (i.e. 1:endd) Fn and Ln data
  #bootstrapping groundphase 
      Fsub<-Fn[1:(dmc[3]-1)]
      Lsub<-Ln[1:(dmc[3]-1)]
      #fitted values are the original fitted values (since the bqseline was optimized to yield residuals with expected value zero in this part)
      #fitg<-Bert(Fsub,ose)
 #preventing problems: remove possible NAs from analysis
      non<-which(is.na(Lsub))
      if(length(non)!=0){
        Lsub<-Lsub[-non]
        Fsub<-Fsub[-non]}
    #Determine cutoff frequency: 
      testfreqs<-round(seq(from=0.1,to=8,length=length(Fsub)*2),digits=2)
      resg<-sapply(testfreqs,lowpass,y=Lsub,RSoS=T)  
      cutfreq<-tail(testfreqs[which(resg>(RSS[1]/4))],1) 
    #apply smoothing (and add fitted values behind the groundphase ones 
      fitg<-lowpass(Lsub,cutfreq,RSoS=F)
    #re-insert NA values
      if(length(non)!=0){for(i in 1:length(non)){
        if(non[i]==1){
          Lsub<-c(NA,Lsub)
          fitg<-c(NA,fitg)
        }else{
          Lsub<-c(Lsub[1:(non[i]-1)],NA,Lsub[non[i]:length(Lsub)])
          fitg<-c(fitg[1:(non[i]-1)],NA,fitg[non[i]:length(fitg)])}}}
    #the X matrix
      Fsub<-Fn[1:(dmc[3]-1)]      
      xix<-matrix(c(rep(1,times=length(Fsub)),Fsub,Fsub^2),nrow=length(Fsub),ncol=3,byrow=F)
    #residuals*scaling HC2 (based on leverage) 
      resg<-(Lsub-fitg)*1/sqrt((1-diag(xix%*%solve(t(xix)%*%xix)%*%t(xix))))
  #boostrapping first phase 
      Fsub<-Fn[dmc[3]:dmc[4]]
      Lsub<-Ln[dmc[3]:dmc[4]]
    #Determine cutoff frequency:
      testfreqs<-round(seq(from=0.1,to=8,length=length(Fsub)*2),digits=2)
      res1<-sapply(testfreqs,lowpass,y=Lsub,RSoS=T)  
      cutfreq<-tail(testfreqs[which(res1>(RSS[2]/4))],1) 
    #apply smoothing (and add fitted values behind the groundphase ones 
      fit1<-lowpass(Lsub,cutfreq,RSoS=F)
    #residuals*scaling HC2 (based on leverage)  (and add residual values behind the groundphase ones 
      res1<-(Lsub-fit1)*1/sqrt((1-diag(mad.hatter(Fsub))))
  #merge both phases
      ress<-append(resg,res1)
      fits<-append(fitg,fit1)}

if(block==2){ #boostrapping second phase
    #Determine cutoff frequency:
      testfreqs<-round(seq(from=0.1,to=8,length=length(Ln)*2),digits=2)
      ress<-sapply(testfreqs,lowpass,y=Fn,RSoS=T)
      cutfreq<-tail(testfreqs[which(ress>(RSS/4))],1)
    #apply smoothing 
      fits<-lowpass(Fn,cutfreq,RSoS=F)
    #residuals*scaling (based on leverage)
      ress<-(Fn-fits)*sqrt(1/(1-diag(mad.hatter(Ln))))}
       
#actual WILD bootstrap
#############################################  
  straps<-matrix(nrow=length(fits),ncol=n,data=rep(fits,times=n),byrow=F)
  stress<-matrix(nrow=length(fits),ncol=n,data=rep(ress,times=n),byrow=F)*matrix(nrow=length(fits),ncol=n,sample(c(-1,1),size=n*length(fits),replace=T))         #phase 1 scaled residuals * Rademacher resamples
  straps<-straps+stress

#actual jackknive
#############################################
  Jn<-length(fits)
  if(block==1){knives<-matrix(nrow=Jn,ncol=Jn,data=rep(Ln[1:dmc[4]],times=Jn),byrow=F)}
  if(block==2){knives<-matrix(nrow=Jn,ncol=Jn,data=rep(Fn,times=Jn),byrow=F)} 
  for(a in 1:Jn){knives[a,a]<-NA}

#for debugging
if(isTRUE(lice)){
      if(block==1){
       x11(14,7);par(mfrow=c(1,2))
       plot(Ln[1:dmc[4]]~Fn[1:dmc[4]],bty="L",main="Bootstrap fit")
        lines(fits~Fn[1:dmc[4]],lty=2,col="orange")
        curve(Bert(x,ose),0,Fn[dmc[4]],lty=2,col="green3",add=T)
       plot(ress,bty="L",main="Bootstrap residuals")
        abline(h=0,col="red")
       dev.set(dev.cur()-1)}  
       if(block==2){  }}


################################################################################
                ########################################
                ##      End of Bootstrap section      ##
                ##      Start of Fitting Section      ##
                ########################################
#################################################################################
                  
#Fit linear model to bootstraps and jackknives + run Kalman Filter
#############################################

  #analyse bootstrap
    if(block==1){Blps<-sapply(c(1:n),function(a){return(ph1fit(straps[dmc[3]:dmc[4],a]))})}

    if(block==2){Bpar2<-sapply(c(1:n),function(a){tlmp<-try(lm(straps[,a]~Ln),silent=TRUE) ;if(inherits(tlmp,"try-error")){
                                                 tlmp<-c(NA,NA)}else{tlmp<-coef(tlmp)      ;if(tlmp[2]>=0){
                                                 B.negative<-function(zz){return(sum((straps[,n]-(zz[1]-zz[2]^2*Ln))^2,na.rm=TRUE))}
                                                 tlmp<-optim(par=c(straps[1,n],0),fn=B.negative)$par; tlmp[2]<--tlmp[2]^2}
                                                 tlmp<-c((-tlmp[1]/tlmp[2]),(1/tlmp[2]))}   ;return(tlmp)})} 

  #analyse jackknive                             #jackknives currently only handles vertical distances:
    if(block==1){Jlps<-sapply(c(1:Jn),function(a){return(ph1fit(knives[dmc[3]:dmc[4],a]))})}
    #ok:                                              
    if(block==2){Jpar2<-sapply(c(1:Jn),function(a){tlmp<-try(lm(knives[,a]~Ln),silent=TRUE) ;if(inherits(tlmp,"try-error")){
                                                  tlmp<-c(NA,NA)}else{tlmp<-coef(tlmp)      ;if(tlmp[2]>=0){
                                                  B.negative<-function(zz){return(sum((knives[,n]-(zz[1]-zz[2]^2*Ln))^2,na.rm=TRUE))}
                                                  tlmp<-optim(par=c(knives[1,n],0),fn=B.negative)$par; tlmp[2]<--tlmp[2]^2}
                                                  tlmp<-c((-tlmp[1]/tlmp[2]),(1/tlmp[2]))}   ;return(tlmp)})} 
#Construct bilinear vectors
#############################################
               
  if(block==2){Bbps<-sapply(c(1:n) ,function(a){return(bpmaker(mp=Bpar1[,a],mp2=Bpar2[,a]))})
               #for the jackknives things are somewhat less straightforward: to cover the sequential omision of all datapoints
               #for the construcion of the biliniear model we should also make the first phase jackknives into bilinear models
               #(using the ose 2nd phase), to that we attach the bilinear models form the second phase jackknives (using the ose 1st phase)
               Jbps<-sapply(c(1:dim(Jpar1)[2]),function(a){return(bpmaker(mp=Jpar1[,a],mp2=ose[4:5]))})
               Jbps<-cbind(Jbps,sapply(c(1:Jn),function(a){return(bpmaker(mp=ose[1:3] ,mp2=Jpar2[,a]))}))}
                 

################################################################################
                ########################################
                ##       End of Fitting Section       ##
                ##       Start of Output Section      ##
                ########################################
################################################################################
  if(block==1){return(list("blp"=Blps,"jlp"=Jlps))}
  if(block==2){return(list("bbp"=Bbps,"jbp"=Jbps))}
#################
}#end of function
         
