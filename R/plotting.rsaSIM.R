#' @title Plots the response surface of a simulated polynomial model
#'
#' @description
#' Plots the response surface of a model based on simulation parameter values.
#'
#' @details
#' Based on polynomial parameters provided by the user, this function simulates data and plots the corresponding response surface, including reversals and accelerations (when applicable). 
#'
#'
#' @name plotting.rsaSIM
#' @param int Intercept (b0)
#' @param x X coefficient (b1)
#' @param y Y coefficient (b2)
#' @param x2 X^2 coefficient (b3)
#' @param xy XY interaction coefficient (b4)
#' @param y2 Y^2 coefficient (b5)
#' @param y3 Y^3 coefficient (b6)
#' @param x2y X^2Y coefficient (b7)
#' @param xy2 XY^2 coefficient (b8)
#' @param x3 X^3 coefficient (b9)
#' @param corr_xy Correlation between predictors (default to 0).
#' @param e_label If "none", no extrema coordinates are projected. Defaults to NULL.
#' @param center Method for centering the predictor variables before the analysis. Default option ("variablewise") centers the predictor variables on \emph{their respective} sample mean. "none" applies no centering. "pooled" centers the predictor variables on their \emph{pooled} sample mean. You should think carefully before applying the "pooled" option, as centering or reducing the predictor variables on common values (e.g., their grand means and SDs) can affect the commensurability of the predictor scales.
#' @param scale Method for scaling the predictor variables before the analysis. Default option ("variablewise") scales the predictor variables on \emph{their respective} sample SD. "none" applies no scaling. "pooled" scales the predictor variables on their \emph{pooled} sample SD. You should think carefully before applying the "pooled" option, as scaling the predictor variables on common values (e.g., their grand SDs) can affect the commensurability of the predictor scales.
#' @param seed Randomization seed for reproducible results
#' @param xlim Limits of the x axis
#' @param ylim Limits of the y axis
#' @param zlim Limits of the z axis
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param zlab Label for z axis
#' @param main the main title of the plot
#' @param cex.main Factor for main title size
#' @param surface Method for the calculation of the surface z values. "predict" takes the predicted values from the model, "smooth" uses a thin plate smoother (function \code{Tps} from the \code{fields} package) of the raw data
#' @param lambda lambda parameter for the smoother. Default (NULL) means that it is estimated by the smoother function. Small lambdas around 1 lead to rugged surfaces, big lambdas to very smooth surfaces.
#' @param rotation Rotation of the 3d surface plot (when type == "3d")
#' @param label.rotation Rotation of the axis labls (when type == "3d")
#' @param gridsize Number of grid nodes in each dimension
#' @param bw Print surface in black and white instead of colors?
#' @param legend Print color legend for z values?
#' @param cex.tickLabel Font size factor for tick labels
#' @param cex.axesLabel Font size factor for axes labels
#' @param type \code{3d} for 3d surface plot, \code{2d} for 2d contour plot, "interactive" for interactive rotatable plot. 
#' @param points A list of parameters which define the appearance of the raw scatter points: 
#'	\itemize{
#'		\item data: Data frame which contains the coordinates of the raw data points. First column = x, second = y, third = z. This data frame is automatically generated when the plot is based on a fitted RSA-object
#'		\item n_random *Number of randomly drawn data points to be plotted from data. Used to avoid cluttering in large datasets.
#'		\item show = TRUE: Should the original data points be overplotted?
#'		\item color = "black": Color of the points. Either a single value for all points, or a vector with the same size as data points provided. If parameter \code{fill} is also defined, \code{color} refers to the border of the points.
#'		\item fill = NULL: Fill of the points. Either a single value for all points, or a vector with the same size as data points provided. As a default, this is set to NULL, which means that all points simply have the color \code{color}.
#' 		\item value="raw": Plot the original z value, "predicted": plot the predicted z value
#'		\item jitter = 0: Amount of jitter for the raw data points. For z values, a value of 0.005 is reasonable
#'		\item cex = .5: multiplication factor for point size. Either a single value for all points, or a vector with the same size as data points provided.
#' 		\item stilt: Should stilts be drawn for selected data points (i.e., lines from raw data points to the floor)? A logical vector with the same size as data points provided, indicating which points should get a stilt.
#' 		\item out.mark = FALSE: If set to TRUE, outliers according to Bollen & Jackman (1980) are printed as red X symbols, but only when they have been removed in the RSA function: \code{RSAmodel(..., out.rm=TRUE)}.
#'			\itemize{
#'				\item If out.rm == TRUE (in RSAmodel()) and out.mark == FALSE (in plotting.rsa()), the outlier is removed from the model and *not plotted* in RSAplot.
#'				\item If out.rm == TRUE (in RSAmodel()) and out.mark == TRUE (in plotting.rsa()), the outlier is removed from the model but plotted and marked in RSAplot.
#'				\item If out.rm == FALSE (in RSAmodel()): Outliers are not removed and cannot be plotted.
#'				\item Example syntax: \code{plotting.rsa(r1, points=list(show=TRUE, out.mark=TRUE))}
#'		}
#'	}
#' As a shortcut, you can also set \code{points=TRUE} to set the defaults.

#' @param model If x is an RSA object: from which model should the response surface be computed?
#' @param acceleration Rates of accelerations along the LOC and LOIC to be inspected (0< |rate| < 1). Passed on internally to \code{ident.ext} 
#' @param FAST If FALSE, will also project response over LOC and LOIC. If TRUE, will only project the response surface (faster option).
#' @param n_sample Size of simulated sample
#' @param demo Do not change that parameter (internal use only)
#' @param fit Do not change that parameter (internal use only)
#' @param param Should the surface parameters a1 to a5 be shown on the plot? In case of a 3d plot a1 to a5 are printed on top of the plot; in case of a contour plot the principal axes are plotted. Surface parameters are not printed for cubic surfaces.
#' @param coefs Should the regression coefficients b1 to b5 (b1 to b9 for cubic models) be shown on the plot? (Only for 3d plot)
#' @param axes *A vector of strings specifying the axes that should be plotted. Can be any combination of c("LOC", "LOIC","r1_LOC","r2_LOC","r1_LOIC","r2_LOIC","a1_LOC","a2_LOC","a1_LOIC","a2_LOIC"). LOC = line of congruence, LOIC = line of incongruence, r1_LOC = first reversal point on LOC,r2_LOC = second reversal point on LOC,r1_LOIC = first reversal point on LOIC,r2_LOIC = second reversal point on LOIC,a1_LOC = first acceleration point on LOC,a2_LOC = second acceleration point on LOC,a1_LOIC = first acceleration point on LOIC,a2_LOIC = second acceleration point on LOIC.
#' @param axesStyles *Define the visual styles of the axes LOC,LOIC,r1_LOC,r2_LOC,r1_LOIC,r2_LOIC,a1_LOC,a2_LOC,a1_LOIC,a2_LOIC. Provide a named list: \code{axesStyles=list(LOC = list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue"))}. It recognizes three parameters: \code{lty}, \code{lwd}, and \code{col}. If you define a style for an axis, you have to provide all three parameters, otherwise a warning will be shown.
#' @param project *A vector of graphic elements that should be projected on the floor of the cube. Can include any combination of c("LOC", "LOIC", "contour", "points"). Note that projected elements are plotted in the order given in the vector (first elements are plotted first and overplotted by later elements).
#' @param maxlines Should the maximum lines be plotted? (red: maximum X for a given Y, blue: maximum Y for a given X). Works only in type="3d"
#' @param link Link function to transform the z axes. Implemented are "identity" (no transformation; default), "probit", and "logit"
#' @param suppress.surface Should the surface be suppressed (only for \code{type="3d"})? Useful for only showing the data points, or for didactic purposes (e.g., first show the cube, then fade in the surface).
#' @param suppress.box Should the surrounding box be suppressed (only for \code{type="3d"})?
#' @param suppress.grid Should the grid lines be suppressed (only for \code{type="3d"})?
#' @param suppress.ticklabels Should the numbers on the axes be suppressed (only for \code{type="3d"})?
#' @param border Should a thicker border around the surface be plotted? Sometimes this border leaves the surrounding box, which does not look good. In this case the border can be suppressed by setting \code{border=FALSE}.
#' @param contour A list defining the appearance of contour lines (aka. height lines). show=TRUE: Should the contour lines be plotted on the 3d wireframe plot? (Parameter only relevant for \code{type="3d"}). color = "grey40": Color of the contour lines. highlight = c(): A vector of heights which should be highlighted (i.e., printed in bold). Be careful: the highlighted line is not necessarily exactly at the specified height; instead the nearest height line is selected.
#' @param hull Plot a bag plot on the surface (This is a bivariate extension of the boxplot. 50\% of points are in the inner bag, 50\% in the outer region). See Rousseeuw, Ruts, & Tukey (1999).
#' @param showSP Plot the stationary point? (only relevant for \code{type="contour"})
#' @param showSP.CI Plot the CI of the stationary point? (only relevant for \code{type="contour"})
#' @param distance A vector of three values defining the distance of labels to the axes
#' @param tck A vector of three values defining the position of labels to the axes (see ?wireframe)
#' @param pal A palette for shading.
#' @param pal.range Should the color range be scaled to the box (\code{pal.range = "box"}, default), or to the min and max of the surface (\code{pal.range = "surface"})? If set to "box", different surface plots can be compared along their color, as long as the zlim is the same for both.
#' @param pad Pad controls the margin around the figure (positive numbers: larger margin, negative numbers: smaller margin)
#' @param claxes.alpha Alpha level that is used to determine the axes K1 and K2 that demarcate the regions of significance for the cubic models "CL" and "RRCL"
#' @param ... Additional parameters passed to the plotting function (e.g., sub="Title"). A useful title might be the R squared of the plotted model: \code{sub = as.expression(bquote(R^2==.(round(getPar(x, "r2", model="CUBIC"), 3))))}
#' @return A plot of the simulated response surface 
#'
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The Bagplot: A Bivariate Boxplot. The American Statistician, 53(4), 382-387. doi:10.1080/00031305.1999.10474494
#'  
#' Núñez-Regueiro, F., Juhel, J. (2022). \emph{Model-Building Strategies in Response Surface Analysis} Manuscript submitted for publication.
#'  
#' Núñez-Regueiro, F., Juhel, J. (2024a). \emph{Response Surface Analysis for the Social Sciences I: Identifying Best-Fitting Polynomial Solutions} Manuscript submitted for publication.
#'  
#' Núñez-Regueiro, F., Juhel, J. (2024b). \emph{Response Surface Analysis for the Social Sciences II: Combinatory Rationales for Complex Polynomial Models} Manuscript submitted for publication.
#' 
#' @seealso \code{\link{plotting.ext}}, \code{\link{RSAmodel}}
#'
#' @examples
#' ######SIMULATE RESPONSE SURFACE OF A MODEL FM26_PARALLELASYMWEAK (b9=1/3*b8)
#' ###### WITH COR(X,Y)=.20
#' ##Define polynomial parameters and correlation
#' b0 <- 0
#' b1 <- -0.064
#' b2 <- -b1*6
#' b3 <- -0.058
#' b5 <- b3*0.7
#' b4 <- 0
#' b6 <- 0
#' b7 <- 0
#' b8 <- -0.10
#' b9 <- 1/3*b8
#' corrXY <- 0.20
#' ##Simulate corresponding surface
#' plotSIM <- plotting.rsaSIM(x=b1,y=b2,x2=b3,xy=b4,y2=b5,x3=b6,x2y=b7,xy2=b8,y3=b9,corr_xy= corrXY)
#' plotSIM 
#' @rdname plotting.rsaSIM
#' @export

plotting.rsaSIM <- function(int=0,x, y,x2,xy,y2, x3, x2y, xy2, y3, 
	corr_xy=0,e_label=NULL,center="variablewise",scale="variablewise",seed= 123, type="3d", model=c("CUBIC","QUADRATIC","FM3_ADDITIVE"),acceleration=c(0,0),FAST=TRUE,n_sample=10000, 
	xlim=c(-3,3), ylim=c(-3,3), zlim=c(-3,3), 
	xlab=NULL, ylab=NULL, zlab=NULL, main="",
	surface="predict", lambda=NULL, 
	suppress.surface=FALSE, suppress.box = FALSE, suppress.grid = FALSE,
	suppress.ticklabels=FALSE,
	rotation=list(x=-63, y=32, z=15), label.rotation=list(x=19, y=-40, z=92), 
	gridsize=21, bw=FALSE, legend=TRUE, param=TRUE, coefs=FALSE,
	axes=c("LOC", "LOIC","r1_LOC","r2_LOC","r1_LOIC","r2_LOIC","a1_LOC","a2_LOC","a1_LOIC","a2_LOIC"),
	axesStyles=list(
		LOC = list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue")),
		LOIC= list(lty="solid",  lwd=2, col=ifelse(bw==TRUE, "black", "blue")),
		PA1 = list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30")),
		PA2 = list(lty="dotted", lwd=2, col=ifelse(bw==TRUE, "black", "gray30")),
		r1_LOC=list(lty="solid",  lwd=2,col=ifelse(bw==TRUE, "black", "green")),
		r2_LOC=list(lty="solid",  lwd=2,col=ifelse(bw==TRUE, "black", "green")),
		r1_LOIC=list(lty="solid",  lwd=2,col=ifelse(bw==TRUE, "black", "red")),
		r2_LOIC=list(lty="solid",  lwd=2,col=ifelse(bw==TRUE, "black", "red")),
		a1_LOC=list(lty="twodash",  lwd=2,col=ifelse(bw==TRUE, "black", "green")),
		a2_LOC=list(lty="twodash",  lwd=2,col=ifelse(bw==TRUE, "black", "green")),
		a1_LOIC=list(lty="twodash",  lwd=2,col=ifelse(bw==TRUE, "black", "red")),
		a2_LOIC=list(lty="twodash",  lwd=2,col=ifelse(bw==TRUE, "black", "red"))

	),
	project=c("contour"), maxlines=FALSE,
	cex.tickLabel=1, cex.axesLabel=1, cex.main=1, 
	points = list(jitter=0.1,show=T, value="predicted"),
	fit=NULL, link="identity", 
	tck=c(1.5, 1.5, 1.5), distance=c(1.3, 1.3, 1.4), border=FALSE, 
	contour = list(show=FALSE, color="grey40", highlight = c()),
	hull=NA, showSP=FALSE, showSP.CI=FALSE, 
	pal=NULL, pal.range="box", 
	pad=0, claxes.alpha=0.05, demo=FALSE, ...) {
	



########SIMULATE PREDICTOR DATA

set.seed(seed)

###Data multivariate normal
mu_x <- mean(xlim)
mu_y <- mean(ylim)
sd_x <- diff(xlim) / 6  # ~99% within limits
sd_y <- diff(ylim) / 6
rho <- 0.20
var_cov <- matrix(c(sd_x^2,rho * sd_x * sd_y,rho * sd_x * sd_y, sd_y^2), nrow = 2)
sim_data <- MASS::mvrnorm(n = n_sample, mu = c(mu_x, mu_y), Sigma = var_cov)
colnames(sim_data) <- c("varx", "vary")
df_simRSA <- as.data.frame(sim_data)

###Predictor standardization
    if (center == "variablewise" ) {
        df_simRSA[, "varx"] <- scale(df_simRSA[, "varx"], center = TRUE, scale = FALSE)
        df_simRSA[, "vary"] <- scale(df_simRSA[, "vary"], center = TRUE, scale = FALSE)
    }
    if (scale == "variablewise" ) {
        df_simRSA[, "varx"] <- scale(df_simRSA[, "varx"], center = FALSE, scale = TRUE)
        df_simRSA[, "vary"] <- scale(df_simRSA[, "vary"], center = FALSE, scale = TRUE)
    }
    if (center == "pooled") {
        pooled.mean <- mean(c(df_simRSA[, "varx"], df_simRSA[, "vary"]))
        df_simRSA[, "varx"] <- df_simRSA[, "varx"] - pooled.mean
        df_simRSA[, "vary"] <- df_simRSA[, "vary"] - pooled.mean
    }
    if (scale == "pooled") {
        pooled.sd <- sd(c(df_simRSA[, "varx"], df_simRSA[, "vary"]))
        df_simRSA[, "varx"] <- df_simRSA[, "varx"]/pooled.sd
        df_simRSA[, "vary"] <- df_simRSA[, "vary"]/pooled.sd
    }

### Generate polynomial terms
df_simRSA$varx2 <- df_simRSA$varx^2
df_simRSA$vary2 <- df_simRSA$vary^2
df_simRSA$varxy <- df_simRSA$varx*df_simRSA$vary
df_simRSA$varx3 <- df_simRSA$varx^3
df_simRSA$vary3 <- df_simRSA$vary^3
df_simRSA$varx2y <- df_simRSA$varx^2*df_simRSA$vary
df_simRSA$varxy2 <- df_simRSA$vary^2*df_simRSA$varx

### Simulate outcome based on polynomial effects
#Sim parameters
b1 <- x
b2 <- y
b3 <- x2
b4 <- xy
b5 <- y2
b6 <- x3
b7 <- x2y
b8 <- xy2
b9 <- y3
b0 <- int

#Deterministic
df_simRSA$varz <- b1 * df_simRSA$varx +
                  b2 * df_simRSA$vary +
                  b3 * df_simRSA$varx2 +
                  b4 * df_simRSA$varxy +
                  b5 * df_simRSA$vary2 +
                  b6 * df_simRSA$varx3 +
                  b7 * df_simRSA$varx2y +
                  b8 * df_simRSA$varxy2 +
                  b9 * df_simRSA$vary3+
                  b0

#Noise (signal to noise ratio= 10)
noise_sd <- sd(df_simRSA$varz)/10
df_simRSA$varz <- df_simRSA$varz+ rnorm(nrow(df_simRSA), 0, noise_sd)



###########ESTIMATE MODEL
#Surface type
surface <- as.character(ifelse(all(c(b6, b7, b8, b9) == 0) & all(c(b3, b4, b5) == 0), "FM3_ADDITIVE",
            ifelse(all(c(b6, b7, b8, b9) == 0), "QUADRATIC", "CUBIC")))


RSA_sim  <- RSAtools::RSAmodel(formula= varz ~ varx*vary, data= df_simRSA, model= surface,scale=F,center=F,estimator="ML",out.rm=F)

###########PLOT SIMULATION MODEL
#Variable names
			if (is.null(xlab)) {
					XLAB <- "varx"
				}
			else{ XLAB <- xlab}	
			if (is.null(ylab)) {
					YLAB <- "vary"
				}
			else{ YLAB <- ylab}	
			if (is.null(zlab)) {
					ZLAB <- "Z"
				}
			else{ ZLAB <- zlab}	

#Surface type
surface <- as.character(ifelse(all(c(b6, b7, b8, b9) == 0) & all(c(b3, b4, b5) == 0), "LINEAR",
            ifelse(all(c(b6, b7, b8, b9) == 0), "QUADRATIC", "CUBIC")))

if(surface=="CUBIC"){
plotRESP <- RSAtools::plotting.rsa(RSA_sim,model= "CUBIC",projet=NULL,type=type,FAST=FAST,n_randomsample=100,acceleration= acceleration,showSP= showSP, showSP.CI= showSP.CI,project=project,points=list(jitter=0.1,show=F, value="predicted"),axesStyles= axesStyles, axes= axes,main= main,coefs = coefs,param = param,legend = legend,distance = c(1.2, 1.2, 1.2),cex.tickLabel = cex.tickLabel, cex.main = cex.main,xlab= XLAB,ylab= YLAB,zlab= ZLAB,cex.axesLabel = cex.axesLabel,xlim= xlim,ylim= ylim,zlim= zlim,hull= hull)
plotEXT <- RSAtools::plotting.ext(RSA_sim,model="CUBIC",names_xLOC=paste(XLAB,"=",YLAB,sep=""),e_label=e_label,names_xLOIC=paste(XLAB,"=-",YLAB,sep=""),names_z=ZLAB,xlim= xlim,zlim= zlim,acceleration= acceleration)
plotEXT[["3D"]] <- plotRESP
#Extrema values
plotEXT[["Values"]] <- RSAtools::ident.ext(RSA_sim,model="CUBIC",acceleration,n_sample=100)
return(plotEXT)
}


if(surface=="QUADRATIC"){
#Plot 
plotRESP <- RSAtools::plotting.rsa(RSA_sim,model= "QUADRATIC",projet=NULL,type=type,FAST=FAST,n_randomsample=100,acceleration= acceleration,showSP= showSP, showSP.CI= showSP.CI,project=project,points=list(jitter=0.1,show=F, value="predicted"),axesStyles= axesStyles, axes= axes,main= main,coefs = coefs,param = param,legend = legend,distance = c(1.2, 1.2, 1.2),cex.tickLabel = cex.tickLabel, cex.main = cex.main,xlab= XLAB,ylab= YLAB,zlab= ZLAB,cex.axesLabel = cex.axesLabel,xlim= xlim,ylim= ylim,zlim= zlim,hull= hull)
plotEXT <- RSAtools::plotting.ext(RSA_sim,model="QUADRATIC",names_xLOC=paste(XLAB,"=",YLAB,sep=""),e_label=e_label,names_xLOIC=paste(XLAB,"=-",YLAB,sep=""),names_z=ZLAB,xlim= xlim,zlim= zlim,acceleration= acceleration)
plotEXT[["3D"]] <- plotRESP
#Extrema values
plotEXT[["Values"]] <- RSAtools::ident.ext(RSA_sim,model="QUADRATIC",acceleration,n_sample=100)
return(plotEXT)
}


if(surface=="LINEAR"){
#Plot 
plotRESP <- RSAtools::plotting.rsa(RSA_sim,model= "FM3_ADDITIVE",projet=NULL,type=type,FAST=FAST,n_randomsample=100,acceleration= acceleration,showSP= showSP, showSP.CI= showSP.CI,project=project,points=list(jitter=0.1,show=F, value="predicted"),axesStyles= axesStyles, axes= axes,main= main,coefs = coefs,param = param,legend = legend,distance = c(1.2, 1.2, 1.2),cex.tickLabel = cex.tickLabel, cex.main = cex.main,xlab= XLAB,ylab= YLAB,zlab= ZLAB,cex.axesLabel = cex.axesLabel,xlim= xlim,ylim= ylim,zlim= zlim,hull= hull)
plotEXT <- RSAtools::plotting.ext(RSA_sim,model="FM3_ADDITIVE",names_xLOC=paste(XLAB,"=",YLAB,sep=""),e_label=e_label,names_xLOIC=paste(XLAB,"=-",YLAB,sep=""),names_z=ZLAB,xlim= xlim,zlim= zlim,acceleration= acceleration)
plotEXT[["3D"]] <- plotRESP
#Extrema values
plotEXT[["Values"]] <- RSAtools::ident.ext(RSA_sim,model="FM3_ADDITIVE",acceleration,n_sample=100)
return(plotEXT)
}



}