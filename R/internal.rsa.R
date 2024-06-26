# internal.rsa.R
# = Internal functions used across various RSAtools functions

#' @keywords internal 
.which.zeroes <- function(x,z){
min_zero <- min(abs(diff(diff(x$z))))
zero <- x[which( abs(diff(diff(x$z))) == min_zero ),]
as.numeric(rownames(zero))
}


#' @keywords internal 
  ####Function which.peaks (for cubic function)
.which.peaks <- function(x,partial=TRUE,decreasing=FALSE){
if (decreasing){
if (partial){
peak <- which(diff(c(TRUE,diff(x)<=0,FALSE))>0)
}else {
peak <- which(diff(diff(x)<=0)>0)
} 
peak[1]
}
else {
if (partial){
peak <- which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
}else {
peak <- which(diff(diff(x)>=0)<0)
}
}
}

#' @keywords internal 
.which.pits <- function(x,partial=TRUE,decreasing=FALSE){
if (decreasing){
if (partial){
which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
}else {
which(diff(diff(x)>=0)<0)
} 
}
else {
if (partial){
which(diff(c(TRUE,diff(x)<=0,FALSE))<0)
}else {
which(diff(diff(x)<=0)<0)
}
}
}


#' @keywords internal 
   .ci_pred2 <- function(obj, x, y, side, n, p, alpha, model){

  DV <- obj$DV
  
  # compute predicted outcome value at (x,y)
  z <- .predictRSA(obj, x, y, model=model)

  # vector of predictor values
  r = c(1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3)
  
  # get covariances of the vector beta of estimated coefficients
  COV0 = vcov(obj$models[[model]])
  covBeta = COV0[c(paste0(DV,"~1"), "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9"),][,c(paste0(DV,"~1"), "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9")]
  
  # Compute standard error of the predicted value z
  var_z = t(r) %*% covBeta %*% r
  se_z = sqrt(var_z)
  
  #NEW: compute two-sided confidence interval of z. 
  if(side=="two-sided"){
    lower_z = as.vector(z - qt(1-alpha/2, df=n-p)*se_z)
    upper_z = as.vector(z + qt(1-alpha/2, df=n-p)*se_z)
  }


  #NOT USED: compute one-sided confidence interval of z. 
  if(side=="right"){
    lower_z = as.vector(z - qt(1-alpha, df=n-p)*se_z)
    # lower_z = as.vector(-Inf)
    upper_z = as.vector(z + qt(1-alpha, df=n-p)*se_z)
    
  }
  
  if(side=="left"){
    lower_z = as.vector(z - qt(1-alpha, df=n-p)*se_z)
    upper_z = as.vector(Inf)
    upper_z = as.vector(z + qt(1-alpha, df=n-p)*se_z)
  }
  
  # output
  out <- list(z=z, lower_z=lower_z, upper_z=upper_z)
  return(out)
  }


#' @keywords internal 
.add.variables <- function(formula, df) {
	IV1 <- all.vars(formula)[2]
	IV2 <- all.vars(formula)[3]	
	IV12 <- paste0(IV1, "2")
	IV22 <- paste0(IV2, "2")
	IV13 <- paste0(IV1, "3")
	IV23 <- paste0(IV2, "3")
	IV_IA <- paste0(IV1, "_", IV2)
	IV_IA2 <- paste0(IV1, "2", "_", IV2)
	IV_IA3 <- paste0(IV1, "_", IV2, "2")		
	df[, IV12] <- df[, IV1]^2
	df[, IV22] <- df[, IV2]^2
	df[, IV_IA] <- df[, IV1]*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - Edwards (2002) model
	df$W.JRE <- ifelse(df[, IV1] >= df[, IV2], 0, 1)
	df[, paste0("W.JRE_", IV1)] <- df$W.JRE*df[, IV1]
	df[, paste0("W.JRE_", IV2)] <- df$W.JRE*df[, IV2]
	
	# three new variables for piecewise regression (test absolute difference score) - new model Schoenbrodt 2012
	df$W <- ifelse(df[, IV1] >= df[, IV2], 1, -1)
	df$W[df[, IV1] == df[, IV2]] <- 0
	df[, paste0("W_", IV1)] <- df$W*df[, IV1]
	df[, paste0("W_", IV2)] <- df$W*df[, IV2]
	
	df$diff <- df[, IV2] - df[, IV1]
	df$SD <- df$diff^2
	df$absdiff <- abs(df$diff)
	
	# cubic terms
	df[, IV13] <- df[, IV1]^3
	df[, IV_IA2] <- df[, IV1]^2*df[, IV2]
	df[, IV_IA3] <- df[, IV1]*df[, IV2]^2
	df[, IV23] <- df[, IV2]^3
	
	return(df)
}


#' @keywords internal 
.anovaList <- function(modellist) {
    mods <- modellist[!sapply(modellist, function(x) is.null(x))]
    mods <- mods[!sapply(mods, function(x) !lavaan::inspect(x, "converged"))]
    if (length(mods) == 0) {
        return(list(n.mods = 0))
    }
    DF <- sapply(mods, lavaan::fitmeasures, "df")
    mods <- mods[order(DF, decreasing = FALSE)]
    mods.scaled <- unlist(lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus", 
            "mean.var.adjusted", "scaled.shifted") %in% unlist(sapply(slot(x, 
            "test"), "[", "test")))
    }))
    if (!(all(mods.scaled) | !any(mods.scaled))) {
        mods[[which(sapply(mods, fitmeasures, "df") == 0)]]@test[[2]]$test <- mods[[which(mods.scaled)[1]]]@test[[2]]$test
    }
    pStr <- sapply(1:length(mods), function(x) {
        if (x == 1) {
            paste("mods[[", x, "]]", sep = "")
        }
        else {
            paste("force(mods[[", x, "]])", sep = "")
        }
    })
    pStr2 <- paste0("lavaan::lavTestLRT(", paste(pStr, collapse = ", "), 
        ", method='default')")
    a1 <- eval(parse(text = pStr2))
    if (length(mods) > 1) {
        rownames(a1) <- names(mods)
    }
    attr(a1, "n.mods") <- length(mods)
    return(list(ANOVA = a1, models = mods, n.mods = length(mods)))
}


#' @keywords internal 
.cModels <- function(mL, set, free.max) {
	aL1 <- .anovaList(mL)
	if (aL1$n.mods > 1) {
		N <- lavaan::nobs(aL1$models[[1]])
		a1 <- cbind(aL1$ANOVA[, c(1, 4:7)], plyr::ldply(aL1$models, function(X) {
			FF <- lavaan::fitmeasures(X)
			R <- lavaan::inspect(X, "r2")
			names(R) <- "R2"
                n <- lavaan::nobs(X)
                k <- free.max - FF["df"]
                suppressWarnings({
                  R2.p <- ifelse(k == 0, NA, pf(((n - k - 1) * 
                    R)/(k * (1 - R)), k, n - k - 1, lower.tail = FALSE))
                })
                names(R2.p) <- "R2.p"
                K <- k + 2
                AICc <- -2 * FF["logl"] + 2 * K + 2 * (K * (K + 
                  1))/(n - K - 1)
                names(AICc) <- NULL
                return(c(AICc = AICc, FF[c("cfi", "srmr")], R, 
                  R2.p))
            }))
        a1 <- a1[, !grepl(".id", colnames(a1))]
        a1$k <- free.max - a1$Df
        a1$R2.adj <- 1 - ((1 - a1$R2)) * ((N - 1)/(N - a1$k - 
            1))
        a1$delta.R2 <- c(NA, a1$R2[1:(nrow(a1) - 1)] - a1$R2[2:(nrow(a1))])
        a1$model <- rownames(a1)
        a1$set <- set
        return(a1)
	}
}



#' @keywords internal 
.R2difftest <- function(x, unrestricted="", restricted="interceptonly"){
  
  n <- lavaan::nobs(x$models[[unrestricted]])
  free.max <- .getFreeParameters(x$models[[unrestricted]])
  
  Fu <- lavaan::fitmeasures(x$models[[unrestricted]])
  Ru <- lavaan::inspect(x$models[[unrestricted]], "r2")
  ku <- free.max - Fu["df"]
  
  if ( restricted=="interceptonly" ){
    Rr <- 0
    kr <- 0
  } else {
    Fr <- lavaan::fitmeasures(x$models[[restricted]])
    Rr <- lavaan::inspect(x$models[[restricted]], "r2")
    kr <- free.max - Fr["df"]
  }
  
  suppressWarnings({		
    R2.p <- ifelse(ku==kr,
                   NA,
                   pf( ( ( n - ku - 1 ) * ( Ru - Rr ) ) / ( (ku-kr) * (1 - Ru ) ), ku-kr, n-ku-1, lower.tail=FALSE)
    )
  })
  
  delta.R2 <- Ru-Rr
  names(delta.R2) <- "R2"
  names(R2.p) <- "R2.p"
  
  return(list(
    delta.R2=delta.R2,
    R2.p=R2.p, 
    ku=ku
  ))
  
}


#' @keywords internal 
.f2 <- function(x, digits=2, prepoint=0, skipZero=FALSE) {
	
	if (skipZero == TRUE) {zero <- "."} else {zero <- "0."}
	
	if (length(dim(x)) == 2) {
		apply(x, 2, function(x2) {gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x2) , fixed=TRUE)})
	} else {
		gsub("0.", zero, sprintf(paste("%",prepoint,".",digits,"f",sep=""), x) , fixed=TRUE)
	}
}

#' @keywords internal 
.p2star <- function(val) {
	
	res <- val
	
	for (i in 1:length(val)) {
		res[i] <- ""
		if (is.na(val[i])) next();
		if (val[i] <= 0.1) res[i] <- "\U2020"
		if (val[i] <= 0.05) res[i] <- "*"
		if (val[i] <= 0.01) res[i] <- "**"
		if (val[i] <= 0.001) res[i] <- "***"
	}
	
	return(res)
}

#' @keywords internal 
.p0 <- function(x) {
	if (is.na(x)) return("NA")
	if (x >= .001) return(paste0("p = ", .f2(x, 3, skipZero=TRUE)))
	if (x <  .001) return("p <.001")	
}
.p <- Vectorize(.p0)

#' @keywords internal 
.getFreeParameters <- function(model) {
	VARS <- nrow(lavaan::inspect(model, "free")$beta)	# number of variables
	df.max <- (VARS*(VARS+1))/2		# maximum df
	df.pred <- ((VARS-1)*(VARS))/2 + 1 # df bound in the predictors (i.e., (co)variances of the predictors & variance of DV)
	free.max <- df.max - df.pred	# maximum of free parameters
	return(free.max)
}


#' @keywords internal 
.getIntersect <- function(b0=0, x=0, y=0, x2=0, xy=0, y2=0, p0, p1, xlim=c(-2, 2), grid=21) {
	X <- seq(min(xlim), max(xlim), length.out=grid)
	Y <- p0 + p1*X
	n <- data.frame(X, Y)
	n2 <- .add.variables(z~X+Y, n)
	n2$Z <- b0 + colSums(c(x, y, x2, y2, xy)*t(n2[, c(1:5)]))
	return(n2[, c("X", "Y", "Z")])
}

#' @keywords internal 
.model <- function(x, model="full") x$models[[model]]
#' @keywords internal 
.syntax <- function(x, model="full") cat(x$models[[model]]@Options$model)



#' @keywords internal 
.pRamp <- function(p, sig=.05, borderline=.10, bias=.8) {
	# calculate bias that the color transition is at the borderline value
	bias2 <- .33/(borderline/(1 - sig))
	cR1 <- colorRamp(c("red", "red", "orange"), bias=bias, space="Lab")
	cR2 <- colorRamp(c("orange", "green", "green"), bias=bias2, space="Lab")
	
	p2 <- rep("#FFFFFF", length(p))
	if (length(p[p < sig])>0) {
		p2[p < sig] <- rgb(cR1(p[p < sig]/sig), maxColorValue=255)
	}
	if (length(p[p >= sig])>0) {
		p2[p >= sig] <- rgb(cR2((p[p >= sig] - sig) / (1 - sig)), maxColorValue=255)
	}
	return(p2)
}


#' @keywords internal 
.f0 <- function (vec, target, unique = TRUE) {
    ret <- vec[sapply(target, function(x) which.min(abs(x - vec)))]
    if (unique) { ret <- unique(ret) }
    ret
}


#' @keywords internal 
.predictRSA <- function(RSAobject, X, Y, model="full") {
	C <- coef(RSAobject$models[[model]])
	if (RSAobject$models[[model]]@Options$estimator != "DWLS") {
		b0 <- as.numeric(ifelse(is.na(C[paste0(RSAobject$DV, "~1")]), b0, C[paste0(RSAobject$DV, "~1")]))
		} else {
			# the threshold is the negative of the intercept ...
			b0 <- -as.numeric(ifelse(is.na(C[paste0(RSAobject$DV, "|t1")]), b0, C[paste0(RSAobject$DV, "|t1")]))
		}
	x <- as.numeric(ifelse(is.na(C["b1"]), 0, C["b1"]))
	y <- as.numeric(ifelse(is.na(C["b2"]), 0, C["b2"]))
	x2 <- as.numeric(ifelse(is.na(C["b3"]), 0, C["b3"]))
	y2 <- as.numeric(ifelse(is.na(C["b5"]), 0, C["b5"]))
	xy <- as.numeric(ifelse(is.na(C["b4"]), 0, C["b4"]))
	w <- as.numeric(ifelse(is.na(C["w1"]), 0, C["w1"]))
	wx <- as.numeric(ifelse(is.na(C["w2"]), 0, C["w2"]))
	wy <- as.numeric(ifelse(is.na(C["w3"]), 0, C["w3"]))
	
	# cubic parameters
	x3 <- as.numeric(ifelse(is.na(C["b6"]), 0, C["b6"]))
	x2y <- as.numeric(ifelse(is.na(C["b7"]), 0, C["b7"]))
	xy2 <- as.numeric(ifelse(is.na(C["b8"]), 0, C["b8"]))
	y3 <- as.numeric(ifelse(is.na(C["b9"]), 0, C["b9"]))
	
	
	C <- c(x, y, x2, y2, xy, w, wx, wy,x3, x2y, xy2, y3)
	
	# compute predicted value
	Z <- b0 + colSums(C*t(cbind(X, Y, X^2, Y^2, X*Y, 0, 0, 0, X^3, X^2*Y, X*Y^2, Y^3)))
	return(Z)
}



#' @keywords internal 
.interpolatePolygon <- function(x, y, minDist, plot=FALSE) {
	minDist <- minDist^2	# compare with squared x^2 + y^2 (faster)
	interp <- data.frame()
	pol <- data.frame(x, y)
	colnames(pol) <- c("x", "y")
	for (i in 1:(nrow(pol)-1)) {
		# get distance
		D <- (pol[i, 1] - pol[i+1, 1])^2 + (pol[i, 2] - pol[i+1, 2])^2
		if (D > minDist) { 
			N <- ceiling(sqrt(D)/sqrt(minDist)) # number of interpolations
			APPROX <- data.frame(
				x = seq(pol[i, 1], pol[i+1, 1], length.out=N),
				y = seq(pol[i, 2], pol[i+1, 2], length.out=N)
			)
			interp <- rbind(interp, APPROX)
		} else if (D>0 & D <= minDist){
			interp <- rbind(interp, pol[i, ])
			if (i==1) colnames(interp) <- c("x", "y")
		}
	}
	interp <- rbind(interp, pol[nrow(pol), ])
	if (plot==TRUE) {
		plot(pol, col="red")
		points(interp[, 1], interp[, 2], col="green", pch=20)
		lines(interp[, 1], interp[, 2], col="darkgreen")
		text(pol, label=1:nrow(pol))
	}
	return(interp)
}



#' @keywords internal 
.compare2 <- function (x, m1 = "", m2 = "full", digits = 3, verbose = TRUE) 
{
    if (is.null(x$models[[m1]]) | is.null(x$models[[m2]])) {
        stop("You need two models for comparison! At least one of the models has not been fit.")
    }
    free.max <- .getFreeParameters(x$models[[m1]])
    mL <- list(M1 = x$models[[m1]], M2 = x$models[[m2]])
    names(mL) <- c(m1, m2)
    res <- .cModels(mL, set = "two_models", free.max)
    if (verbose == TRUE & !is.null(res)) {
        print(round(res[, 1:13], digits))
    }
    invisible(res)
}


