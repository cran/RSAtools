#' @title Probe extrema in the response surface
#'
#' @description
#' Identify reversal or acceleration points (generically called "extrema") in the LOC or LOIC of the response surface and test how many of them have outcome observations that significantly differ from what would be expected for predictor combinations on these points (that have the same level)
# '
#' @details
#' When testing for reversals or accelerations in nonlinear response surfaces involving quadratic or cubic polynomial families (FM4 to FM37), the \code{RSAextrema} function helps to determine the exact location of reversal or acceleration points along the lines of congruenc (LOC) or incongruence (LOIC), and the number and percentage of observations significantly affected by these reversals or accelerations (for a given probability level, alpha). This points are determined according to derivatives of the function according to rationales for combining polynomials (Núñez-Regueiro & Juhel, 2022, 2024).
#' 
# '
# '
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats qt
# '
#' @export
#' @param RSA_object An RSA object
#' @param model The model to be probed for extrema (reversal or acceleration points)
#' @param acceleration Rates of accelerations along the LOC and LOIC to be inspected (0< abs(rate) < 1). Acceleration points will only appear if reversals do not exist, and if acceleration rates exist (if not, a warning will appear). 
#' @param alpha Alpha level for the one-sided confidence interval of the outcome predictions on the extrema
#' @param z_tested Should significance tests be conducted on "observed" or "predicted" observation
#' @param alphacorrection Set "Bonferroni" to adjust the alpha level for multiple testing when testing the outcome predictions of all data points behind the extrema
#' @param n_sample Number of random draws to consider to find extrema. This option is used for large samples to increase speed in preliminary analyses, but it is not recommended for published results). Defaults to NULL.
#' @param verbose Should extra information be printed?
#' @param df_out Number of random draws to consider to find extrema. This option is used for large samples to increase speed in preliminary analyses, but it is not recommended for published results). Defaults to NULL.
#' 
#' @return A table containing the location and percentages of observations above or below extrema
#' 
#' @references 
#' Núñez-Regueiro, F., Juhel, J. (2022). \emph{Model-Building Strategies in Response Surface Analysis} Manuscript submitted for publication.
#'  
#' Núñez-Regueiro, F., Juhel, J. (2024). \emph{Response Surface Analysis for the Social Sciences II: Combinatory Rationales for Complex Polynomial Models} Manuscript submitted for publication.
#' 
#' @seealso \code{\link{plotting.ext}}, \code{\link{RSAmodel}}
#'
#' @export

ident.ext <- function(RSA_object,model = NULL,acceleration=c(0,0), alpha = 0.05, z_tested="observed", alphacorrection = "none",n_sample=NULL,verbose = TRUE,  df_out=FALSE) {

### Function
    rsa = RSA_object
    if (is.null(model))
    	stop("Please specify which model from the RSA_object should be estimated") 
    if (lavaan::inspect(rsa$models[[model]], "converged") == FALSE) {
        warning("The model has not converged!")
        return(NULL)
    }
    DV <- rsa$DV
    IV1 <- rsa$IV1
    IV2 <- rsa$IV2
    C <- coef(rsa$models[[model]])
    b0 = as.numeric(coef(rsa$models[[model]])[paste0(DV, "~1")])
	c1 = as.numeric(coef(rsa$models[[model]])["b3"])
    c2 = as.numeric(coef(rsa$models[[model]])["b6"])

####Adjust intercept for control variables (if any)
 	if(rsa$is.cv){
	  # vector of means of the control variables
	  cvmeans <- colMeans(rsa$data[, rsa$control.variables], na.rm=T)
	  # coefficients of the control variables
	  cvcoefs <- C[paste0(rsa$DV, "~", rsa$control.variables)]
	  # new b0 = b0 + cv1*mean(controlvariable1) + cv2*mean(controlvariable2) + ...
	  b0 <- b0 + sum(cvmeans*cvcoefs)	  
	  }
 
####Select data
    data = rsa$data
    if(!is.null(n_sample)){
    	# set.seed(123)
   	data <- data[sample(nrow(data), n_sample),] 	
    }
    
    df = data[data$out == FALSE, ]


####Select outcome (observed, predicted) to be compared with reversal/acceleration points
if(z_tested=="observed" | is.null(z_tested) ){
	df$z.ind <- df[,DV]
}
if(z_tested=="predicted"){
    df$z.ind <- NA
    for (k in 1:nrow(df)) {
        df[k, "z.ind"] = .predictRSA(rsa, df[k, IV1], df[k, IV2], 
            model = model)
    }
 }


 #### Transform cluster var to numeric to avoid error in .ci_pred2
# if(length(rsa$models[[model]]@Data@cluster)!=0  & class(rsa$models[[model]]@Data@cluster)!="factor"){
cluster_var <- rsa$models[[model]]@Data@cluster		
if(length(cluster_var)!=0  & inherits(cluster_var,"factor")==F ){
name_cluster <- rsa$models[[model]]@Data@cluster
var_cluster <- df[, name_cluster]
df[, name_cluster] <- as.numeric(as.factor(df[, name_cluster]))
}


 
 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### Compute extreme points on the LOC and the LOIC

 #Output matrix
row_extrema <- c("r1_LOC","r2_LOC","r1_LOIC","r2_LOIC","a1_LOC","a2_LOC","a1_LOIC","a2_LOIC")
col_extrema <- c("X_value","Y_value","Z_value","Pct_Above","Pct_AboveandDiffZ","Pct_Below","Pct_BelowandDiffZ")
matrix_ext <- matrix(nrow=length(row_extrema),ncol=length(col_extrema),dimnames=list(row_extrema, col_extrema))

 #Data range and model
x_min <-  min(data[,c(IV1)],na.rm=T)
x_max <-  max(data[,c(IV1)],na.rm=T)
xlim <- c(x_min, x_max)
var_x <- seq(xlim[1],xlim[2], psych::describe(seq(xlim[1],xlim[2]))$se/5000)
C <- coef(rsa$models[[model]])
param_model <- lavaan::parameterEstimates(rsa$models[[model]])

 #### #### #### LOC  #### #### ####
 #Coef LOC
u1 <- ( C["b1"]+ C["b2"] )
u2 <- ( C["b3"]+C["b4"]+C["b5"] )
u3 <- (C["b6"]+C["b7"]+C["b8"]+C["b9"])

####Correct non-zero estimates of true zeros (to avoid invalid values in roots)
u1 <- ifelse(abs(u3/u1)>10^4,0,u1)	 
u2 <- ifelse(u2!=0 & abs(u3/u2)>10^4,0,u2) 	
u3 <- ifelse(u3!=0 &abs(u1/u3)>10^4,0,u3) 		
u3 <- ifelse(u2!=0 & abs(u2/u3)>10^4,0,u3) 

####Roots for z_LOC (when applicable)
if(u3!=0){	
x_r1_LOC =  -( ((u2^2)-3*u1*u3)^0.5 +u2)/(3*u3)
x_r2_LOC =  ( ((u2^2)-3*u1*u3)^0.5 -u2)/(3*u3)
}

if(u3==0 | round(u3,8)==0 ){
x_r1_LOC =  -u1/(2*u2)
x_r2_LOC =  NA
}

#NA impossible roots (applies when u3 is zero but estimated !=0)
impossible_E <-  abs(max(c(xlim))+3*sd(var_x))
x_r1_LOC <- ifelse(abs(x_r1_LOC) > impossible_E,NA, x_r1_LOC)
x_r2_LOC <- ifelse(abs(x_r2_LOC) > impossible_E,NA, x_r2_LOC)

#Response value on roots
z_r1_LOC <- b0 +u1*x_r1_LOC + u2*x_r1_LOC ^2 + u3*x_r1_LOC ^3 
z_r2_LOC <- b0 +u1*x_r2_LOC + u2*x_r2_LOC ^2 + u3*x_r2_LOC ^3 

#Approximation of roots and respone values on roots using simulated data (for case where u3 is fixed to zero but estimated at != 0, the calculation of roots is flawed)
z_LOC <- b0 +u1*var_x + u2* var_x^2 + u3* var_x^3 
df_LOC <- data.frame(var_x, z_LOC)
peaks_LOC <- .which.peaks(df_LOC$z_LOC, decreasing=FALSE)
r1_LOC <- ifelse(u3 < 0 & length(peaks_LOC) > 1, peaks_LOC[2],peaks_LOC[1] )
r1_LOC <- ifelse(df_LOC[c(r1_LOC),"var_x"] == max(var_x)   |  df_LOC[c(r1_LOC),"var_x"] == min(var_x) | abs(df_LOC[c(r1_LOC),"var_x"]) > impossible_E, NA, r1_LOC )
pits_LOC <- .which.pits(df_LOC$z_LOC, decreasing=FALSE)
r2_LOC <- ifelse(u3 < 0 | length(pits_LOC) == 1, pits_LOC[1], pits_LOC[2] )
r2_LOC <- ifelse(df_LOC[c(r2_LOC),"var_x"] == max(var_x)  | df_LOC[c(r2_LOC),"var_x"] == min(var_x) | abs(df_LOC[c(r2_LOC),"var_x"]) > impossible_E, NA, r2_LOC )
extremes_LOC <- df_LOC[c(r1_LOC, r2_LOC),]
extremes_LOC <- na.omit(extremes_LOC)

# Reorder e1 and e2 from small to high x
if(!is.na(extremes_LOC[1,"var_x"]) & !is.na(extremes_LOC[2,"var_x"]) ){
if(extremes_LOC[1,"var_x"] > extremes_LOC[2,"var_x"]){
xr1_LOC <- extremes_LOC[1,]	
xr2_LOC <- extremes_LOC[2,]	
extremes_LOC[1,]	<- xr2_LOC
extremes_LOC[2,]	<- xr1_LOC
}
}
if( !is.na(x_r1_LOC) & !is.na(x_r2_LOC)){
if(x_r1_LOC > x_r2_LOC){
ir1_LOC <- r1_LOC
ir2_LOC <- r2_LOC
r1_LOC <- ir2_LOC
r2_LOC <- ir1_LOC
xr1_LOC <- x_r1_LOC	
xr2_LOC <- x_r2_LOC	
x_r1_LOC	<- xr2_LOC
x_r2_LOC	<- xr1_LOC
zr1_LOC <- z_r1_LOC	
zr2_LOC <- z_r2_LOC	
z_r1_LOC	<- zr2_LOC
z_r2_LOC	<- zr1_LOC
}
}


### Extrema on LOC : x = -y
if(!is.na(x_r1_LOC) ){
matrix_ext["r1_LOC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(x_r1_LOC,3),round(x_r1_LOC,3),round(z_r1_LOC,3))))
}
if( is.na(param_model[which(param_model$label=="u3"),"se"]) && is.na(x_r1_LOC) && !is.na(r1_LOC) && round(u3,8)==0 && round(u2,8)!=0){
matrix_ext["r1_LOC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOC[c(r1_LOC),"var_x"],3),round(df_LOC[c(r1_LOC),"var_x"],3),round(df_LOC[c(r1_LOC),"z_LOC"],3))))
}
if(!is.na(x_r2_LOC)){
matrix_ext["r2_LOC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(x_r2_LOC,3),round(x_r2_LOC,3),round(z_r2_LOC,3))))
}
if( is.na(param_model[which(param_model$label=="u3"),"se"]) && is.na(x_r2_LOC) && !is.na(r2_LOC) && round(u3,8)==0 && round(u2,8)!=0){
matrix_ext["r2_LOC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOC[c(r2_LOC),"var_x"],3),round(df_LOC[c(r2_LOC),"var_x"],3),round(df_LOC[c(r2_LOC),"z_LOC"],3))))
}




 #### #### #### LOIC  #### #### ####
 #Coef LOIC
v1 <- ( C["b1"]- C["b2"])
v2 <- ( C["b3"]-C["b4"]+C["b5"] )
v3 <- (C["b6"]-C["b7"]+C["b8"]-C["b9"])

####Correct non-zero estimates of true zeros (to avoid invalid values in roots)
v1 <- ifelse(abs(v3/v1)>10^4,0,v1)	 
v2 <- ifelse(v2!=0 & abs(v3/v2)>10^4,0,v2) 	
v3 <- ifelse(v3!=0 &abs(v1/v3)>10^4,0,v3) 		
v3 <- ifelse(v2!=0 & abs(v2/v3)>10^4,0,v3) 

####Roots for z_LOIC (when applicable)
if(v3!=0){	
x_r1_LOIC =  -( ((v2^2)-3*v1*v3)^0.5 +v2)/(3*v3)
x_r2_LOIC =  ( ((v2^2)-3*v1*v3)^0.5 -v2)/(3*v3)
}

if(v3==0 | round(v3,8)==0 ){
x_r1_LOIC =  -v1/(2*v2)
x_r2_LOIC =  NA
}


#NA impossible roots (applies when v3 is zero but estimated !=0)
impossible_E <-  abs(max(c(xlim))+3*sd(var_x))
x_r1_LOIC <- ifelse(abs(x_r1_LOIC) > impossible_E,NA, x_r1_LOIC)
x_r2_LOIC <- ifelse(abs(x_r2_LOIC) > impossible_E,NA, x_r2_LOIC)

#Response value on roots
z_r1_LOIC <- b0 +v1*x_r1_LOIC + v2*x_r1_LOIC ^2 + v3*x_r1_LOIC ^3 
z_r2_LOIC <- b0 +v1*x_r2_LOIC + v2*x_r2_LOIC ^2 + v3*x_r2_LOIC ^3 

#Approximation of roots and response values on roots using simulated data (for cases where v3 is fixed to zero but estimated at != 0, the computation of roots is flawed)
z_LOIC <- b0 +v1*var_x + v2* var_x^2 + v3* var_x^3 
df_LOIC <- data.frame(var_x, z_LOIC)
peaks_LOIC <- .which.peaks(df_LOIC$z_LOIC, decreasing=FALSE)
r1_LOIC <- ifelse(v3 < 0 & length(peaks_LOIC) > 1, peaks_LOIC[2],peaks_LOIC[1] )
r1_LOIC <- ifelse(df_LOIC[c(r1_LOIC),"var_x"] == max(var_x)   |  df_LOIC[c(r1_LOIC),"var_x"] == min(var_x) | abs(df_LOIC[c(r1_LOIC),"var_x"]) > impossible_E, NA, r1_LOIC )
pits_LOIC <- .which.pits(df_LOIC$z_LOIC, decreasing=FALSE)
r2_LOIC <- ifelse(v3 < 0 | length(pits_LOIC) == 1, pits_LOIC[1], pits_LOIC[2] )
r2_LOIC <- ifelse(df_LOIC[c(r2_LOIC),"var_x"] == max(var_x)  | df_LOIC[c(r2_LOIC),"var_x"] == min(var_x) | abs(df_LOIC[c(r2_LOIC),"var_x"]) > impossible_E, NA, r2_LOIC )
extremes_LOIC <- df_LOIC[c(r1_LOIC, r2_LOIC),]
extremes_LOIC <- na.omit(extremes_LOIC)

# Reorder e1 and e2 from small to high x
if(!is.na(extremes_LOIC[1,"var_x"]) & !is.na(extremes_LOIC[2,"var_x"]) ){
if(extremes_LOIC[1,"var_x"] > extremes_LOIC[2,"var_x"]){
xr1_LOIC <- extremes_LOIC[1,]	
xr2_LOIC <- extremes_LOIC[2,]	
extremes_LOIC[1,]	<- xr2_LOIC
extremes_LOIC[2,]	<- xr1_LOIC
}
}
if( !is.na(x_r1_LOIC) & !is.na(x_r2_LOIC)){
if(x_r1_LOIC > x_r2_LOIC){
ir1_LOIC <- r1_LOIC
ir2_LOIC <- r2_LOIC
r1_LOIC <- ir2_LOIC
r2_LOIC <- ir1_LOIC
xr1_LOIC <- x_r1_LOIC	
xr2_LOIC <- x_r2_LOIC	
x_r1_LOIC	<- xr2_LOIC
x_r2_LOIC	<- xr1_LOIC
zr1_LOIC <- z_r1_LOIC	
zr2_LOIC <- z_r2_LOIC	
z_r1_LOIC	<- zr2_LOIC
z_r2_LOIC	<- zr1_LOIC
}
}

### Extrema on LOIC : x = -y
if(!is.na(x_r1_LOIC) ){
matrix_ext["r1_LOIC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(x_r1_LOIC,3),-round(x_r1_LOIC,3),round(z_r1_LOIC,3))))
}
if(is.na(param_model[which(param_model$label=="v3"),"se"]) && is.na(x_r1_LOIC) && !is.na(r1_LOIC) && round(v3,8)==0 && round(v2,8)!=0){
matrix_ext["r1_LOIC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOIC[c(r1_LOIC),"var_x"],3),-round(df_LOIC[c(r1_LOIC),"var_x"],3),round(df_LOIC[c(r1_LOIC),"z_LOIC"],3))))
}
if(!is.na(x_r2_LOIC) ){
matrix_ext["r2_LOIC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(x_r2_LOIC,3),-round(x_r2_LOIC,3),round(z_r2_LOIC,3))))
}
if(is.na(param_model[which(param_model$label=="v3"),"se"]) && is.na(x_r2_LOIC) && !is.na(r2_LOIC) && round(v3,8)==0 && round(v2,8)!=0){
matrix_ext["r2_LOIC",c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOIC[c(r2_LOIC),"var_x"],3),-round(df_LOIC[c(r2_LOIC),"var_x"],3),round(df_LOIC[c(r2_LOIC),"z_LOIC"],3))))
}


matrix_ext <- data.frame(matrix_ext)


 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### If E1 does not exist but E2 does, E2= E1
if( (is.na(matrix_ext["r1_LOC",c("X_value")]) && !is.na(matrix_ext["r2_LOC",c("X_value")])) |(!is.na(matrix_ext["r2_LOC",c("X_value")]) & matrix_ext["r1_LOC",c("X_value")]==matrix_ext["r2_LOC",c("X_value")]) ){
matrix_ext["r1_LOC",] <- matrix_ext["r2_LOC",]
matrix_ext["r2_LOC",] <- NA
}

if( (is.na(matrix_ext["r1_LOIC",c("X_value")]) && !is.na(matrix_ext["r2_LOIC",c("X_value")])) |(!is.na(matrix_ext["r2_LOIC",c("X_value")]) & matrix_ext["r1_LOIC",c("X_value")]==matrix_ext["r2_LOIC",c("X_value")]) ){
matrix_ext["r1_LOIC",] <- matrix_ext["r2_LOIC",]
matrix_ext["r2_LOIC",] <- NA
}



 

 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### Corrections made for close to zero but non-zero extrema (rounded to zero)
 #### NB:This happens when metaparameters are fixed to zero 
 ####	but are estimated by lavaan as non-zero (10^-5)

#LOC
if( round(u3,5)!=0 && round(u2,5)==0 && u3/u1<0   && is.na(matrix_ext[c("r1_LOC","r2_LOC"),c("X_value")]) ){
r1_LOC <-  .which.zeroes(x=df_LOC,z=z_LOC)
r2_LOC <-  .which.zeroes(x=df_LOC,z=z_LOC)
matrix_ext[c("r1_LOC"),c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOC[c(r2_LOC),"var_x"],3),-round(df_LOC[c(r2_LOC),"var_x"],3),round(df_LOC[c(r2_LOC),"z_LOC"],3))))	
 matrix_ext[c("r2_LOC"),c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOC[c(r2_LOC),"var_x"],3),-round(df_LOC[c(r2_LOC),"var_x"],3),round(df_LOC[c(r2_LOC),"z_LOC"],3))))	
 }
 
#LOIC
if( round(v3,5)!=0 && round(v2,5)==0 && v3/v1<0 && is.na(matrix_ext[c("r1_LOIC","r2_LOIC"),c("X_value")]) ){
r1_LOIC <-  .which.zeroes(x=df_LOIC,z=z_LOIC)
r2_LOIC <-  .which.zeroes(x=df_LOIC,z=z_LOIC)
matrix_ext[c("r1_LOIC"),c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOIC[c(r2_LOIC),"var_x"],3),-round(df_LOIC[c(r2_LOIC),"var_x"],3),round(df_LOIC[c(r2_LOIC),"z_LOIC"],3))))	
 matrix_ext[c("r2_LOIC"),c("X_value","Y_value","Z_value")] <- as.numeric(t(c(round(df_LOIC[c(r2_LOIC),"var_x"],3),-round(df_LOIC[c(r2_LOIC),"var_x"],3),round(df_LOIC[c(r2_LOIC),"z_LOIC"],3))))	
 } 
 



 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### Compute acceleration points on the LOC and the LOIC 
 #### for acceleration rates in acceleration[1] and acceleration[2] 

if(round(u3,5)!=0 && round(u1,5)!=0 && u3/u1> 0   && is.na(matrix_ext[c("r1_LOC"),c("X_value")]) && acceleration[1]!=0 ){
### Compute acceleration points (when they exist, i.e., u1/u3 > 0)
#Lowest rate of change
df_LOC$rate_zLOC <- b0 +1*u1+2*var_x*u2+3*var_x^2*u3
x_lowrate_LOC <- df_LOC[which(df_LOC$rate_zLOC== min(df_LOC$rate_zLOC)),"var_x"]; x_lowrate_LOC
rate0_ZLOC <- min(df_LOC$rate_zLOC); rate0_ZLOC

###Rate of change at acceleration factor (default= 2)
rate0_ZLOC_acceleration <- acceleration[1]

###X values at accelerated rate of change
xLOC_e1_acceleration =  -( ( u2^2 -3*(u1-rate0_ZLOC_acceleration)*u3)^0.5 +u2)/(3*u3)
zLOC_e1_acceleration <- b0 +xLOC_e1_acceleration* u1+xLOC_e1_acceleration^2*u2+ xLOC_e1_acceleration^3*u3
xLOC_e2_acceleration =  ( ( u2^2-3*(u1-rate0_ZLOC_acceleration)*u3)^0.5 -u2)/(3*u3)
zLOC_e2_acceleration <- b0 +xLOC_e2_acceleration* u1+xLOC_e2_acceleration^2*u2+ xLOC_e2_acceleration^3*u3
matrix_ext[c("a1_LOC"),c("X_value","Y_value","Z_value")] <- t(cbind(round(xLOC_e1_acceleration,3),round(xLOC_e1_acceleration,3),round(zLOC_e1_acceleration,3)))	
matrix_ext[c("a2_LOC"),c("X_value","Y_value","Z_value")] <- t(cbind(round(xLOC_e2_acceleration,3),round(xLOC_e2_acceleration,3),round(zLOC_e2_acceleration,3)))	

if( is.na(xLOC_e1_acceleration)  ){
stop("The specified acceleration rate along the LOC does not exist. Please reverse the sign or increase the size of the acceleration.")}
if( xLOC_e1_acceleration > xLOC_e2_acceleration ){
SAVEe1 <- matrix_ext[c("a1_LOC"),c("X_value","Y_value","Z_value")] 
SAVEe2 <- matrix_ext[c("a2_LOC"),c("X_value","Y_value","Z_value")] 
matrix_ext[c("a1_LOC"),c("X_value","Y_value","Z_value")] <- SAVEe2
matrix_ext[c("a2_LOC"),c("X_value","Y_value","Z_value")] <- SAVEe1
}

 }

 
if(round(v3,5)!=0 && round(v1,5)!=0 && v3/v1> 0   && is.na(matrix_ext[c("r1_LOIC"),c("X_value")])  && acceleration[2]!=0 ){
### Compute acceleration points (when they exist, i.e., u1/u3 > 0)
#Lowest rate of change
df_LOIC$rate_zLOIC <- b0 +1*v1+2*var_x*v2+3*var_x^2*v3
x_lowrate_LOIC <- df_LOIC[which(df_LOIC$rate_zLOIC== min(df_LOIC$rate_zLOIC)),"var_x"]; x_lowrate_LOIC
rate0_ZLOIC <- min(df_LOIC$rate_zLOIC); rate0_ZLOIC

###Rate of change at acceleration factor (default= 2)
rate0_ZLOIC_acceleration <- acceleration[2]

###X values at accelerated rate of change
xLOIC_e1_acceleration =  -( ( v2^2-3*(v1-rate0_ZLOIC_acceleration)*v3)^0.5 +v2)/(3*v3)
zLOIC_e1_acceleration <- b0 +xLOIC_e1_acceleration*v1+xLOIC_e1_acceleration^2*v2+ xLOIC_e1_acceleration ^3*v3
xLOIC_e2_acceleration =  ( ( v2^2-3*(v1-rate0_ZLOIC_acceleration)*v3)^0.5 -v2)/(3*v3)
zLOIC_e2_acceleration <- b0 +xLOIC_e2_acceleration*v1+xLOIC_e2_acceleration^2*v2+ xLOIC_e2_acceleration ^3*v3
matrix_ext[c("a1_LOIC"),c("X_value","Y_value","Z_value")] <- t(cbind(round(xLOIC_e1_acceleration,3),-round(xLOIC_e1_acceleration,3),round(zLOIC_e1_acceleration,3)))	
matrix_ext[c("a2_LOIC"),c("X_value","Y_value","Z_value")] <- t(cbind(round(xLOIC_e2_acceleration,3),-round(xLOIC_e2_acceleration,3),round(zLOIC_e2_acceleration,3)))	

if( is.na(xLOIC_e1_acceleration)  ){
stop("The specified acceleration rate along the LOIC does not exist. Please reverse the sign or increase the size of the acceleration.")}
if( xLOIC_e1_acceleration > xLOIC_e2_acceleration ){
SAVEe1 <- matrix_ext[c("a1_LOIC"),c("X_value","Y_value","Z_value")] 
SAVEe2 <- matrix_ext[c("a2_LOIC"),c("X_value","Y_value","Z_value")] 
matrix_ext[c("a1_LOIC"),c("X_value","Y_value","Z_value")] <- SAVEe2
matrix_ext[c("a2_LOIC"),c("X_value","Y_value","Z_value")] <- SAVEe1
}

 } 





 
 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### Compute percentage of data points above/below extrema and with signifcantly different response (at alpha level, one-tailed)
 
###NB: The confidence interval to which the predicted value of each data point behind extrema must be compared is specific to the respective data point. (cf. ca_range RRCA, Humberg et al., 2020). Specific reponse values on the line (xi,yi) that intersects with extrema on the LOIC and LOC are defined, respectively, by yi(LOC) = 1*(xi)-(xExt) and yi(LOIC) = -1*(xi)+(xExt)
###NB2: .ci_pred2 works with two-sided alpha	==> e.g., for one-tailed .05, requires .10
alpha2side <- alpha*2

#### r1_LOC CUBIC
if(!is.na(matrix_ext["r1_LOC","X_value"]) & !is.na(matrix_ext["r2_LOC","X_value"]) & round(u3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- -1*df[,IV1]+matrix_ext["r1_LOC","X_value"]

#pct above and below point
df$Below_r1_LOC <- df[,c(IV1)] < df$e_specific_y
# df$Below_r1_LOC[is.na(df$Below_r1_LOC)] <- FALSE
df$Above_r1_LOC <- df[,c(IV1)] > df$e_specific_y
# df$Above_r1_LOC[is.na(df$Above_r1_LOC)] <- FALSE
matrix_ext["r1_LOC","Pct_Below"] <- round(sum(df$Below_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
matrix_ext["r1_LOC","Pct_Above"] <- round(sum(df$Above_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Below_r1_LOC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$BelowAndDiff_r1_LOC <- ifelse((df$Below_r1_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_r1_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["r1_LOC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndDiff_r1_LOC <- ifelse((df$Above_r1_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Above_r1_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["r1_LOC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndLower_r1_LOC <- ifelse( (df$Below_r1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_r1_LOC <- ifelse(  (df$Below_r1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$AboveAndLower_r1_LOC <- ifelse( (df$Above_r1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_r1_LOC <- ifelse(  (df$Above_r1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

} 

#### r2_LOC CUBIC
if(!is.na(matrix_ext["r1_LOC","X_value"]) & !is.na(matrix_ext["r2_LOC","X_value"]) & round(u3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- -1*df[,IV1]+matrix_ext["r2_LOC","X_value"]

#pct above and below point
df$Above_r2_LOC <- df[,c(IV1)] > df$e_specific_y	
# df$Above_r2_LOC[is.na(df$Above_r2_LOC)] <- FALSE
matrix_ext["r2_LOC","Pct_Above"] <- round(sum(df$Above_r2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_r2_LOC <- df[,c(IV1)] < df$e_specific_y	
# df$Below_r2_LOC[is.na(df$Below_r2_LOC)] <- FALSE
matrix_ext["r2_LOC","Pct_Below"] <- round(sum(df$Below_r2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_r2_LOC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})

df$AboveAndDiff_r2_LOC <- ifelse((df$Above_r2_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Above_r2_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["r2_LOC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndDiff_r2_LOC <- ifelse((df$Below_r2_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_r2_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["r2_LOC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndLower_r2_LOC <- ifelse( (df$Above_r2_LOC == T & df$z.ind < df[,"lower_zE"]) , TRUE, FALSE)
df$AboveAndHigher_r2_LOC <- ifelse(   (df$Above_r2_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$BelowAndLower_r2_LOC <- ifelse( (df$Below_r2_LOC == T & df$z.ind < df[,"lower_zE"]) , TRUE, FALSE)
df$BelowAndHigher_r2_LOC <- ifelse(   (df$Below_r2_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

} 


#### r1_LOC QUADRATIC
if(!is.na(matrix_ext["r1_LOC","X_value"]) & (is.na(matrix_ext["r2_LOC","X_value"]) | round(u3,8)==0) ){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- -1*df[,IV1]+matrix_ext["r1_LOC","X_value"]

#pct above and below point
df$Above_r1_LOC <- df[,c(IV1)] > df$e_specific_y
# df$Above_r1_LOC[is.na(df$Above_r1_LOC)] <- FALSE
matrix_ext["r1_LOC","Pct_Above"] <- round(sum(df$Above_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_r1_LOC <- df[,c(IV1)] < df$e_specific_y
# df$Below_r1_LOC[is.na(df$Below_r1_LOC)] <- FALSE
matrix_ext["r1_LOC","Pct_Below"] <- round(sum(df$Below_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_r1_LOC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$AboveAndDiff_r1_LOC <- ifelse( (df$Above_r1_LOC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_r1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndDiff_r1_LOC <- ifelse( (df$Below_r1_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_r1_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)

df$BelowAndLower_r1_LOC <- ifelse( (df$Below_r1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_r1_LOC <- ifelse(  (df$Below_r1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$AboveAndLower_r1_LOC <- ifelse( (df$Above_r1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_r1_LOC <- ifelse(  (df$Above_r1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

matrix_ext["r1_LOC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
matrix_ext["r1_LOC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 


#### r1_LOIC CUBIC
if(!is.na(matrix_ext["r1_LOIC","X_value"]) & !is.na(matrix_ext["r2_LOIC","X_value"]) & round(v3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- df[,IV1]-matrix_ext["r1_LOIC","X_value"]

#pct above and below point
df$Below_r1_LOIC <- df[,c(IV2)] > df$e_specific_y 
# df$Below_r1_LOIC[is.na(df$Below_r1_LOIC)] <- FALSE
matrix_ext["r1_LOIC","Pct_Below"] <- round(sum(df$Below_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Above_r1_LOIC <- df[,c(IV2)] < df$e_specific_y 
# df$Above_r1_LOIC[is.na(df$Above_r1_LOIC)] <- FALSE
matrix_ext["r1_LOIC","Pct_Above"] <- round(sum(df$Above_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Below_r1_LOIC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$BelowAndDiff_r1_LOIC <- ifelse( (df$Below_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Below_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndLower_r1_LOIC <- ifelse( (df$Below_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_r1_LOIC <- ifelse(  (df$Below_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$AboveAndLower_r1_LOIC <- ifelse( (df$Above_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_r1_LOIC <- ifelse(  (df$Above_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

matrix_ext["r1_LOIC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndDiff_r1_LOIC <- ifelse( (df$Above_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
matrix_ext["r1_LOIC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 


#### r2_LOIC CUBIC
if(!is.na(matrix_ext["r1_LOIC","X_value"]) & !is.na(matrix_ext["r2_LOIC","X_value"]) & round(v3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- df[,IV1]-matrix_ext["r2_LOIC","X_value"]

#pct above and below point
df$Above_r2_LOIC <- df[,c(IV2)] < df$e_specific_y 
# df$Above_r2_LOIC[is.na(df$Above_r2_LOIC)] <- FALSE
matrix_ext["r2_LOIC","Pct_Above"] <- round(sum(df$Above_r2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_r2_LOIC <- df[,c(IV2)] > df$e_specific_y 
# df$Below_r2_LOIC[is.na(df$Below_r2_LOIC)] <- FALSE
matrix_ext["r2_LOIC","Pct_Below"] <- round(sum(df$Below_r2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_r2_LOIC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$AboveAndDiff_r2_LOIC <- ifelse( (df$Above_r2_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_r2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndLower_r2_LOIC <- ifelse( (df$Above_r2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_r2_LOIC <- ifelse( (df$Above_r2_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$BelowAndLower_r2_LOIC <- ifelse( (df$Below_r2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_r2_LOIC <- ifelse( (df$Below_r2_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
matrix_ext["r2_LOIC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndDiff_r2_LOIC <- ifelse( (df$Below_r2_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Below_r2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
matrix_ext["r2_LOIC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 


#### r1_LOIC QUADRATIC
if(!is.na(matrix_ext["r1_LOIC","X_value"]) & (is.na(matrix_ext["r2_LOIC","X_value"]) | round(v3,8)==0) ){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- df[,IV1]-matrix_ext["r1_LOIC","X_value"]

#pct above and below point
df$Above_r1_LOIC <- df[,c(IV2)] < df$e_specific_y 
# df$Above_r1_LOIC[is.na(df$Above_r1_LOIC)] <- FALSE
matrix_ext["r1_LOIC","Pct_Above"] <- round(sum(df$Above_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_r1_LOIC <- df[,c(IV2)] > df$e_specific_y 
# df$Below_r1_LOIC[is.na(df$Below_r1_LOIC)] <- FALSE
matrix_ext["r1_LOIC","Pct_Below"] <- round(sum(df$Below_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_r1_LOIC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = "two-sided", n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$AboveAndDiff_r1_LOIC <- ifelse( (df$Above_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndDiff_r1_LOIC <- ifelse( (df$Below_r1_LOIC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_r1_LOIC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
df$BelowAndLower_r1_LOIC <- ifelse( (df$Below_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_r1_LOIC <- ifelse( (df$Below_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

df$AboveAndLower_r1_LOIC <- ifelse( (df$Above_r1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_r1_LOIC <- ifelse( (df$Above_r1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

matrix_ext["r1_LOIC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
matrix_ext["r1_LOIC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_r1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 

matrix_ext$N_total <- rsa[[1]][[model]]@Data@nobs[[1]]





#### #### #### #### #### #### #### #### #### #### #### #### 
 #### Compute percentage of data points above/below acceleration points and with signifcantly different response (at alpha level)
 
#### a1_LOC CUBIC
if(!is.na(matrix_ext["a1_LOC","X_value"]) & !is.na(matrix_ext["a2_LOC","X_value"]) & round(u3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- -1*df[,IV1]+matrix_ext["a1_LOC","X_value"]

#pct above and below point
df$Below_a1_LOC <- df[,c(IV1)] < df$e_specific_y
# df$Below_a1_LOC[is.na(df$Below_a1_LOC)] <- FALSE
df$Above_a1_LOC <- df[,c(IV1)] > df$e_specific_y
# df$Above_a1_LOC[is.na(df$Above_a1_LOC)] <- FALSE
matrix_ext["a1_LOC","Pct_Below"] <- round(sum(df$Below_a1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
matrix_ext["a1_LOC","Pct_Above"] <- round(sum(df$Above_a1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Below_a1_LOC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$BelowAndDiff_a1_LOC <- ifelse((df$Below_a1_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_a1_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["a1_LOC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_a1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndDiff_a1_LOC <- ifelse((df$Above_a1_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Above_a1_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["a1_LOC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_a1_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndLower_a1_LOC <- ifelse( (df$Below_a1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_a1_LOC <- ifelse(  (df$Below_a1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$AboveAndLower_a1_LOC <- ifelse( (df$Above_a1_LOC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_a1_LOC <- ifelse(  (df$Above_a1_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

} 

#### a2_LOC CUBIC
if(!is.na(matrix_ext["a1_LOC","X_value"]) & !is.na(matrix_ext["a2_LOC","X_value"]) & round(u3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- -1*df[,IV1]+matrix_ext["a2_LOC","X_value"]

#pct above and below point
df$Above_a2_LOC <- df[,c(IV1)] > df$e_specific_y	
# df$Above_a2_LOC[is.na(df$Above_a2_LOC)] <- FALSE
matrix_ext["a2_LOC","Pct_Above"] <- round(sum(df$Above_a2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_a2_LOC <- df[,c(IV1)] < df$e_specific_y	
# df$Below_a2_LOC[is.na(df$Below_a2_LOC)] <- FALSE
matrix_ext["a2_LOC","Pct_Below"] <- round(sum(df$Below_a2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_a2_LOC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(u3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})

df$AboveAndDiff_a2_LOC <- ifelse((df$Above_a2_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Above_a2_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["a2_LOC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_a2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndDiff_a2_LOC <- ifelse((df$Below_a2_LOC == T & df$z.ind < df[,"lower_zE"]) | (df$Below_a2_LOC == T & df$z.ind > df[,"upper_zE"]), TRUE, FALSE)
matrix_ext["a2_LOC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_a2_LOC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndLower_a2_LOC <- ifelse( (df$Above_a2_LOC == T & df$z.ind < df[,"lower_zE"]) , TRUE, FALSE)
df$AboveAndHigher_a2_LOC <- ifelse(   (df$Above_a2_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$BelowAndLower_a2_LOC <- ifelse( (df$Below_a2_LOC == T & df$z.ind < df[,"lower_zE"]) , TRUE, FALSE)
df$BelowAndHigher_a2_LOC <- ifelse(   (df$Below_a2_LOC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

} 


#### a1_LOIC CUBIC
if(!is.na(matrix_ext["a1_LOIC","X_value"]) & !is.na(matrix_ext["a2_LOIC","X_value"]) & round(v3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- df[,IV1]-matrix_ext["a1_LOIC","X_value"]

#pct above and below point
df$Below_a1_LOIC <- df[,c(IV2)] > df$e_specific_y 
# df$Below_a1_LOIC[is.na(df$Below_a1_LOIC)] <- FALSE
matrix_ext["a1_LOIC","Pct_Below"] <- round(sum(df$Below_a1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Above_a1_LOIC <- df[,c(IV2)] < df$e_specific_y 
# df$Above_a1_LOIC[is.na(df$Above_a1_LOIC)] <- FALSE
matrix_ext["a1_LOIC","Pct_Above"] <- round(sum(df$Above_a1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Below_a1_LOIC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$BelowAndDiff_a1_LOIC <- ifelse( (df$Below_a1_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Below_a1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndLower_a1_LOIC <- ifelse( (df$Below_a1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_a1_LOIC <- ifelse(  (df$Below_a1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$AboveAndLower_a1_LOIC <- ifelse( (df$Above_a1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_a1_LOIC <- ifelse(  (df$Above_a1_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)

matrix_ext["a1_LOIC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_a1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$AboveAndDiff_a1_LOIC <- ifelse( (df$Above_a1_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_a1_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
matrix_ext["a1_LOIC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_a1_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 


#### a2_LOIC CUBIC
if(!is.na(matrix_ext["a1_LOIC","X_value"]) & !is.na(matrix_ext["a2_LOIC","X_value"]) & round(v3,8)!=0){

#extrema specific to individual
df$e_specific_x <- df[,IV1]
df$e_specific_y <- df[,IV1]-matrix_ext["a2_LOIC","X_value"]

#pct above and below point
df$Above_a2_LOIC <- df[,c(IV2)] < df$e_specific_y 
# df$Above_a2_LOIC[is.na(df$Above_a2_LOIC)] <- FALSE
matrix_ext["a2_LOIC","Pct_Above"] <- round(sum(df$Above_a2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$Below_a2_LOIC <- df[,c(IV2)] > df$e_specific_y 
# df$Below_a2_LOIC[is.na(df$Below_a2_LOIC)] <- FALSE
matrix_ext["a2_LOIC","Pct_Below"] <- round(sum(df$Below_a2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)

#pct above and below point and different Z at 1-alpha CI
alpha_X <- ifelse(alphacorrection == "Bonferroni",alpha2side/sum(df$Above_a2_LOIC,na.rm=T),alpha2side)    
df[,"zE_ci"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$z})
df[,"lower_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$lower_z})
df[,"upper_zE"] <- apply(df, 1, function(data){.ci_pred2(obj = rsa, x = data["e_specific_x"], y = data["e_specific_y"], side = ifelse(v3 < 0, "right", "left"), n = sum(!is.na(df[,c(IV1)]),na.rm=T), p = 10, alpha = alpha_X,model = model)$upper_z})
df$AboveAndDiff_a2_LOIC <- ifelse( (df$Above_a2_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Above_a2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndLower_a2_LOIC <- ifelse( (df$Above_a2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$AboveAndHigher_a2_LOIC <- ifelse( (df$Above_a2_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
df$BelowAndLower_a2_LOIC <- ifelse( (df$Below_a2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
df$BelowAndHigher_a2_LOIC <- ifelse( (df$Below_a2_LOIC == T & df$z.ind > df[,"upper_zE"]) , TRUE, FALSE)
matrix_ext["a2_LOIC","Pct_AboveandDiffZ"] <- round(sum(df$AboveAndDiff_a2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
df$BelowAndDiff_a2_LOIC <- ifelse( (df$Below_a2_LOIC == T & df$z.ind > df[,"upper_zE"]) | (df$Below_a2_LOIC == T & df$z.ind < df[,"lower_zE"]), TRUE, FALSE)
matrix_ext["a2_LOIC","Pct_BelowandDiffZ"] <- round(sum(df$BelowAndDiff_a2_LOIC,na.rm=T)/sum(!is.na(df[,c(IV1)]),na.rm=T)*100,2)
} 





 #### #### #### #### #### #### #### #### #### #### #### #### 
 #### Infer which kind of congruence processes and rationales are at play (for future releases)


# ###CORRECT REDUNDANCIES
if(!is.na(matrix_ext["r2_LOC","X_value"]) & matrix_ext["r1_LOC","X_value"]==matrix_ext["r2_LOC","X_value"]){
matrix_ext["r2_LOC",] <- NA	
}
if(!is.na(matrix_ext["r2_LOIC","X_value"]) & matrix_ext["r1_LOIC","X_value"]==matrix_ext["r2_LOIC","X_value"]){
matrix_ext["r2_LOIC",] <- NA	
}


########################################
#####Return matrix
#Reinject cluster var into df
if(length(cluster_var)!=0  & inherits(cluster_var,"factor")==F ){
df[,name_cluster] <- var_cluster
}
#Get dataframe containing test values (optional)
if(df_out==T){
matrix_ext <- list(EXT=data.frame(matrix_ext),DF_EXT=df)
}
matrix_ext
}

