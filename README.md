
#########################
##  How to use         ##
## Full Process PCR    ##
#########################

Description
	determine the maximal (initial) efficiency of a PCR reaction

Usage
	analyse(data,plots=FALSE,param=FALSE)

Arguments

	data	= a vector containing the consecutive fluorescence measurments.
	plots	= Logical, if TRUE a number of informative graphical representations of the data are returned in addition to the other output
		  the plots are: - data and fitted model (actual fluorescence readings + values reconstructed by the model in green)
				 - Efficiency plots: (A) double logarithmic transformed efficiency vs. fluorescence + bilinear model (vertical dotted lines represent the portion of the data between which the first phase was fitted) 
						     (B) efficiency vs. cycle + FPK model
	param	= Logical, if TRUE the function returns the individual model parameters rather than the standard output
	silent  = Logical, if TRUE no messages are returned 
	vecto	= Logical, if TRUE the data are returned as a vector rather than a list. 


Output
	Standard output returns three values:
		Cq: quantification cycle, calculated at the first positive second derivative maximum of a five parameter logistic model.
		Emax: initial (maximal) efficiency of the reaction
		a*i0: initial reaction fluorescence (number of initial target copies times the fluorescence of a single amplicon), expressed in Fluorescence Units (FU)

	parameter output returns three sets of parameters:

		Baseline: parameters of the baseline (linear model)
			intercept = y-intercept of the linear model
			slope	  = slope of the linear model 

		5PLM: parameters of the five paramter logistic model (with sloped baseline)
			fmax + s*x + (y0 - fmax) / (1 + (2^(1/g)-1) * exp(b*(x-xmid)))^g
			
			y0	= a numeric parameter representing the intercept of the baseline
			s 	= a numeric parameter representing the slope of the baseline
			fmax	= a numeric parameter representing the horizontal asymptote on the right side (plateau)
			xmid	= a numeric parameter representing the input value at the center of the curve
			b	= a numeric scale parameter on the input axis, the ’growth rate’
			g	= the reciprocal of the shape parameter: affects near which asymptote maximum ’growth’ occurs

		Bilinear: parameters of the bilinear model
			 chi + eta *log(exp(( a1a * (input-ipt)^2+ a1b *(input-ipt))/eta) + exp( a2 *(n-ipt)/eta))

			 a1a	= "slope" describes together with a1b the curve of the first phase of efficiency decline
		         a1b	= "slope" describes together with a1a the curve of the first phase of efficiency decline
		          a2	= "slope" represents the steepness of the second phase of efficiency decline
		         ipt	= corresponds to the horizontal position of the phase change
		         eta	= affects the abruptness of transition between the two phases
		         chi	= affects the y-intercept of the curve
				

Example
	react<-c(193.6422, 209.2176,  221.9136,  233.8279,  243.1534,  245.9008,  251.7054,  253.5101,  260.7267,  258.2353,  261.5051,  265.9691,  268.3236,  273.0342,  274.9183, 281.7958,  291.9236,  312.8073,  353.2469,  425.0797,  558.4097,  801.8603, 1141.0918, 1617.8774, 2179.0478, 2796.2294, 3426.8560, 4029.0126, 4564.6015, 5065.4347, 5478.7222, 5865.2510, 6197.9947, 6502.7750, 6783.9353, 6997.8949, 7211.0589, 7422.2636, 7626.5014, 7732.2221, 7820.8930, 7925.0876, 7991.5905, 8029.5963, 8045.4035, 8045.6049, 8064.2315, 8103.9244, 8120.5715, 8104.7643, 8117.3816, 8090.2008, 8132.5172, 8078.8213, 8081.2783, 8058.4114, 8021.9082, 8020.0389, 7986.3723, 7999.9835)
	analyse(react)
	
	#or
	analyse(react,plots=TRUE)

	#or
	analyse(react,param=TRUE)


##############################################################################################################
Disclaimer 
	This function has only been tested on a limited amount of data, there may be a few bugs left. 
	Please contact the author with any remarks and/or suggestions.

	email: antoon.lievens@wiv-isp.be

	Most problems are caused by a lack of data resulting in a suboptimal fit try running more cycles first! 
	60 cycles will do in most cases. Increasing primer concentration may help to increase the length of 
	the firs phase of efficiency decline resulting in a better fit. 

	If problems perist, feel free to contact me.
##############################################################################################################	
