#' @title Estimate polynomial models for response surface analysis
#' @description Estimates 37 polynomial families ("STEP1") or user-specific polynomial models ("USER") for use in response surface analysis, based on a 3-step identification strategy (Núñez-Regueiro & Juhel, 2022, 2024).
#' @details This function implements a comparative framework for identifying best-fitting RSA solutions (Núñez-Regueiro & Juhel, 2022, 2024). The default feature ("STEP1") involves the comparison of 37 polynomial families predefined by parametric constraints, to identify likely candidates for best-fitting solution. Step 2 involves probing variants within the retained best-fitting family, by testing user-specific constraints ("USER") on lower-order polynomials that do not define the family. Step 3 (optional) can be conducted on the final variant by bootstrapping validation (using \code{bootstrapLavaan(RSA_object$models$name_final, FUN="coef")}, see examples) or cross-validation data.
#' @param formula A formula in the form \code{z ~ x*y}, specifying the variable names used from the data frame, where z is the name of the response variable, and x and y are the names of the predictor variables.
#' @param data A data frame with the variables
#' @param center Method for centering the predictor variables before the analysis. Default option ("none") applies no centering. "pooled" centers the predictor variables on their \emph{pooled} sample mean. "variablewise" centers the predictor variables on \emph{their respective} sample mean. You should think carefully before applying the "variablewise" option, as centering or reducing the predictor variables on common values (e.g., their grand means and SDs) can affect the commensurability of the predictor scales.
#' @param scale Method for scaling the predictor variables before the analysis. Default option ("none") applies no scaling. "pooled" scales the predictor variables on their \emph{pooled} sample SD, which preserves the commensurability of the predictor scales. "variablewise" scales the predictor variables on \emph{their respective} sample SD. You should think carefully before applying the "variablewise" option, as scaling the predictor variables at different values (e.g., their respective SDs) can affect the commensurability of the predictor scales.
#' @param na.rm Remove missings before proceeding?
#' @param add Additional syntax that is added to the lavaan model. Can contain, for example, additional constraints, like "p01 == 0; p11 == 0"
#' @param out.rm Should outliers according to Bollen & Jackman (1980) criteria be excluded from the analyses? In large data sets this analysis is the speed bottleneck. If you are sure that no outliers exist, set this option to FALSE for speed improvements.
#' @param breakline Should the breakline in the unconstrained absolute difference model be allowed (the breakline is possible from the model formulation, but empirically rather unrealistic ...). Defaults to \code{FALSE}
#' @param verbose Should additional information during the computation process be printed?
#' @param models A vector with names of all models that should be computed. Should be any from \code{c("CUBIC","FM1_ONLYX","FM2_ONLYY","FM3_ADDITIVE","FM4_INTER","FM5_QUADX","FM6_QUADY","FM7_CONG","FM8_INCONG","FM9_CURVCONGX","FM10_CURVCONGY","FM11_CURVINCONGX","FM12_CURVINCONGY","FM13_QUADXQUADY","FM14_ROTCONG","FM15_ROTINCONG","FM16_CUBICX","FM17_CUBICY","FM18_LEVDEPQUADX","FM19_LEVDEPQUADY","FM20_ASYMCONG","FM21_ASYMINCONG","FM22_LEVDEPCONG","FM23_LEVDEPINCONG","FM24_PARALLELASYM","FM25_NONPARALLELASYM","FM26_PARALLELASYMWEAK","FM27_PARALLELASYMSTRONG","FM28_NONPARALLELASYMWEAK","FM29_NONPARALLELASYMSTRONG","FM30_ASYMCONGROTY","FM31_ASYMCONGROTX","FM32_ASYMINCONGROTY","FM33_ASYMINCONGROTX","FM34_LEVDEPCONGROTY","FM35_LEVDEPCONGROTX","FM36_LEVDEPINCONGROTY","FM37_LEVDEPINCONGROTX")}. For \code{models="STEP1"}, all polynomial families and the saturated cubic are computed (default), for \code{models="USER"} all user-specific models defined in the list "user_model" are computed.
#' @param user_model A list of user-specified polynomial models, defined by setting constraints on the polynomial parameters b1 to b9, using lavaan syntax (e.g., "b1 == 2*b2"). Only parametric constraints are allowed.
#' @param estimator Type of estimator that should be used by lavaan. Defaults to "MLR", which provides robust standard errors, a robust scaled test statistic, and can handle missing values. If you want to reproduce standard OLS estimates, use \code{estimator="ML"} and \code{se="standard"}
#' @param se Type of standard errors. This parameter gets passed through to the \code{sem} function of the \code{lavaan} package. See options there. By default, robust SEs are computed. If you use \code{se="boot"}, \code{lavaan} provides CIs and p-values based on the bootstrapped standard error. If you use \code{confint(..., method="boot")}, in contrast, you get CIs and p-values based on percentile bootstrap.
#' @param missing Handling of missing values (this parameter is passed to the \code{lavaan} \code{sem} function). By default (\code{missing=NA}), Full Information Maximum Likelihood (FIML) is employed in case of missing values. If cases with missing values should be excluded, use \code{missing = "listwise"}.
#' @param control.variables A string vector with variable names from \code{data}. These variables are added as linear predictors to the model (in order "to control for them"). No interactions with the other variables are modeled.
#' @param center.control.variables Should the control variables be centered before analyses? This can improve interpretability of the intercept, which will then reflect the predicted outcome value at the point (X,Y)=(0,0) when all control variables take their respective \emph{average} values.
#' @param sampling.weights Name of variable containing sampling weights. Needs to be added here (not in ...) to be included in the analysis dataset.
#' @param group_name Name of variable defining groups, for multigroup modeling.
#' @param cluster Name of variable for clusters, for cluster-level variance correction.
#' @param ... Additional parameters passed to the \code{\link[lavaan]{sem} } function.
#' @return A list of objects containing polynomials models and names of variables 
#' @references
#' Núñez-Regueiro, F., Juhel, J. (2022). \emph{Model-Building Strategies in Response Surface Analysis.} Manuscript submitted for publication.
#'  
#' Núñez-Regueiro, F., Juhel, J. (2024). \emph{Response Surface Analysis for the Social Sciences I: Identifying Best-Fitting Polynomial Solutions.} Manuscript submitted for publication.
#' 
#' @seealso \code{\link{plotting.rsa}}, \code{\link{plotting.ext}}
#' 
#' @examples
#' ###Simulation data from RSAtools
#' summary(sim_NSfit)
#' 
#' ##### STEP 1: Estimate and compare polynomial families
#' ##NB: to estimate all families, use "STEP1" in the "model" option
#' RSA_step1 <-  RSAmodel(engagement ~ needs*supplies,
#' data= sim_NSfit, model= c("CUBIC","FM20_ASYMCONG",
#' "FM21_ASYMINCONG","FM26_PARALLELASYMWEAK"))
#' names(RSA_step1$models)
#' ##Compare solutions according to parsimony (wAIC=Akaike weight)
#' RSA_step1_fit <- best.rsa(RSA_step1,order=c("wAIC"))
#' head(RSA_step1_fit)
#' ##Inspect best-fitting family model
#' summary(RSA_step1$models$FM26_PARALLELASYMWEAK)
#' 
#' ##### STEP 2: Test variant within best-fitting family
#' ##Define variant as constraints
#' list_variant <- list()
#' list_variant[["variant1"]] <- c('
#' ####First-order polynomials: variant-specific			
#' 			b1 == -1/4*b2
#' ####Second-order polynomials: variant-specific			
#' 			b5 == b3/2
#' 			b4 == 0
#' ####Third-order polynomials: define family FM26			
#' 			b6 == 0
#' 			b7 == 0
#' 			b9 == b8/-3		
#' ')
#' RSA_NSfit_Step2  <- RSAmodel(formula= engagement ~ needs*supplies,
#' data= sim_NSfit, model= c("FM26_PARALLELASYMWEAK","USER"),
#' user_model= list_variant)
#' ##Compare variant to best-fitting family (e.g., LRT_pvalue p > .05)
#' best.rsa2(RSA_NSfit_Step2,m1="variant1",m2="FM26_PARALLELASYMWEAK")[2,1:3]
#' 
#' ######PROBE RESPONSE SURFACE EXTREMA: see ident.ext()
#' 
#' ######PLOT RESPONSE SURFACE: see plotting.rsa, plotting.ext
#' 
#' @export


RSAmodel <- function(formula, data = NULL, center = "none", scale = "none",
    na.rm = FALSE, out.rm = TRUE, breakline = FALSE, models = c("CUBIC","STEP1"),user_model=NULL,
    verbose = TRUE, add = "", estimator = "MLR", se = "robust", missing = NA,
    control.variables = NULL, center.control.variables = FALSE, sampling.weights =NULL,group_name=NULL,cluster=NULL,
    ...)
{


if (!requireNamespace("semTools", quietly = TRUE)) {
			stop('`semTools` package needed for smooth surfaces to work. Please install it with install.packages("semTools")', call. = FALSE)
		}

### is.cubic function
	    cubicmodels <- c("cubic", "CA", "CL", "RRCA", "RRCL")
      is.cubic <- any(models %in% cubicmodels)

### empty model fit
    s.cubic <- s.ONLYX <- s.ONLYY <- s.ADDITIVE <- s.INTER <- s.QUADX <- s.QUADY <- s.CONG <- s.INCONG <- s.CURVCONGX <- s.CURVCONGY <- s.CURVINCONGX <- s.CURVINCONGY <- s.QUADXQUADY <- s.ROTCONG <- s.ROTINCONG <- s.CUBICX <- s.CUBICY <- s.LEVDEPQUADX <- s.LEVDEPQUADY <- s.ASYMCONG <- s.ASYMINCONG <- s.LEVDEPCONG <- s.LEVDEPINCONG <- s.PARALLELASYM <- s.NONPARALLELASYM <- s.PARALLELASYMWEAK <- s.PARALLELASYMSTRONG <- s.NONPARALLELASYMWEAK <- s.NONPARALLELASYMSTRONG <- s.ASYMCONGROTY <- s.ASYMCONGROTX <- s.ASYMINCONGROTY <- s.ASYMINCONGROTX <- s.LEVDEPCONGROTY <- s.LEVDEPCONGROTX <- s.LEVDEPINCONGROTY <- s.LEVDEPINCONGROTX <- s.USER <- s.quadratic  <- NULL



s.USER <- list(USER= s.USER)

    SRSQD.rot <- ""
    SRRR.rot <- ""
    add <- paste0("\n# User defined syntax:\n", add)
    DV <- all.vars(formula)[1]
    IV1 <- all.vars(formula)[2]
    IV2 <- all.vars(formula)[3]
    df <- data[, c(DV, IV1, IV2, control.variables)]
#####NEW: multigroup option
    if (!is.null(group_name)) {
    df <- data[, c(DV, IV1, IV2, control.variables,group_name)]
    }
    if (center == TRUE) {
        center <- "variablewise"
    }
    if (scale == TRUE) {
        scale <- "variablewise"
    }
    if (center == "variablewise") {
        warning("You specified 'variablewise' centering of the predictor variables. Make sure to check whether this is a good choice, as it can distort commensurability. Use center='pooled' if you want to center at the pooled mean instead.",
            call. = FALSE)
        df[, IV1] <- scale(df[, IV1], center = TRUE, scale = FALSE)
        df[, IV2] <- scale(df[, IV2], center = TRUE, scale = FALSE)
    }
    if (scale == "variablewise") {
        warning("You specified 'variablewise' scaling of the predictor variables. Make sure to check whether this is a good choice, as it can distort commensurability. Use scale='pooled' if you want to scale at the pooled SD instead.",
            call. = FALSE)
        df[, IV1] <- scale(df[, IV1], center = FALSE, scale = TRUE)
        df[, IV2] <- scale(df[, IV2], center = FALSE, scale = TRUE)
    }
    if (center == "pooled") {
        pooled.mean <- mean(c(df[, IV1], df[, IV2]), na.rm = T)
        df[, IV1] <- df[, IV1] - pooled.mean
        df[, IV2] <- df[, IV2] - pooled.mean
    }
    if (scale == "pooled") {
        pooled.sd <- sd(c(df[, IV1], df[, IV2]), na.rm = T)
        df[, IV1] <- df[, IV1]/pooled.sd
        df[, IV2] <- df[, IV2]/pooled.sd
    }
    df <- .add.variables(formula, data.frame(data.matrix(df)))
    if (0 < min(df[, IV1], na.rm = TRUE) | 0 > max(df[, IV1],
        na.rm = TRUE))
        warning(paste("The numerical zero point is outside of the range of variable",
            IV1, ". Please consider re-centering the variable."))
    if (0 < min(df[, IV2], na.rm = TRUE) | 0 > max(df[, IV2],
        na.rm = TRUE))
        warning(paste("The numerical zero point is outside of the range of variable",
            IV2, ". Please consider re-centering the variable."))
    if ((max(df[, IV1], na.rm = TRUE) - min(df[, IV1], na.rm = TRUE))/(max(df[,
        IV2], na.rm = TRUE) - min(df[, IV2], na.rm = TRUE)) >
        2)
        warning("Predictor variables have a very different range (by factor 2 or larger) - please check scaling of variables.")
    if (is.na(missing)) {
        if (any(is.na(df))) {
            missing <- "fiml"
            warning("There are missing values in your data set. Model is computed with option `missing = 'fiml'`. This is only valid if the data are missing completely at random (MCAR) or missing at random (MAR)! If you want to exclude NAs, use `missing = 'listwise'`",
                call. = FALSE)
        }
        else {
            missing <- "listwise"
        }
    }
    if (any(is.na(df)) & missing == "fiml" & packageVersion("lavaan") <
        "0.6.3") {
        stop("Please install the latest version of 'lavaan'! lavaan versions < 0.6.3 will apply listwise deletion (because fixed.x = TRUE in the sem() models) and have no workaround implemented, so they should not be applied.")
    }
    if (any(is.na(df)) & missing == "fiml" & packageVersion("lavaan") >=
        "0.6.3") {
        missing <- "fiml.x"
    }
    IV12 <- paste0(IV1, "2")
    IV22 <- paste0(IV2, "2")
    IV13 <- paste0(IV1, "3")
    IV23 <- paste0(IV2, "3")
    IV_IA <- paste0(IV1, "_", IV2)
    IV_IA2 <- paste0(IV1, "2", "_", IV2)
    IV_IA3 <- paste0(IV1, "_", IV2, "2")
    W_IV1 <- paste0("W_", IV1)
    W_IV2 <- paste0("W_", IV2)
    is.cv <- length(control.variables) > 0
    if (is.cv & center.control.variables) {
        df[, control.variables] <- scale(df[, control.variables],
            center = TRUE, scale = FALSE)
    }
    CV <- ifelse(is.cv, paste0(" + ", paste(control.variables,
        collapse = " + ")), "")
    addcubic <- ""
    if (is.cubic)
        addcubic <- paste0(" + ", paste(IV13, IV_IA2, IV_IA3,
            IV23, sep = " + "))
    f <- paste0(paste0(DV, " ~ ", paste(IV1, IV2, IV12, IV_IA,
        IV22, sep = " + ")), addcubic, CV)
    lm.full <- lm(f, df, na.action = na.exclude)
    if (is.null(out.rm) || (typeof(out.rm) == "logical" && out.rm ==
        TRUE)) {
        out.rm <- "bj1980"
    }
    if ((typeof(out.rm) == "logical" && out.rm == FALSE)) {
        out.rm <- "none"
    }
    out.rm <- match.arg(out.rm, c("bj1980", "robust", "none"))
    df$out <- FALSE
    if (out.rm == "bj1980") {
        inf <- influence.measures(lm.full)
        df$out <- apply(inf$is.inf[, c("dffit", "cook.d", "hat")],
            1, sum) == 3
        n.out <- sum(na.omit(df$out) == TRUE)
        if (verbose == TRUE & n.out > 0) {
            warning(paste("Removed", n.out, "multivariate outlier(s) according to Bollen & Jackman (1980) criteria. Outliers are in row(s):",
                paste(which(df$out == TRUE), collapse = ", ")))
        }
    }
    if (out.rm == "robust") {
        stop("Robust outlier detection not implemented yet.")
    }
    df$out[is.na(df$out)] <- FALSE
    withCallingHandlers({
        poly <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2, " + b3*",
            IV12, " + b4*", IV_IA, " + b5*", IV22, CV)
        polycubic <- paste0(DV, " ~ b1*", IV1, " + b2*", IV2,
            " + b3*", IV12, " + b4*", IV_IA, " + b5*", IV22,
            " + b6*", IV13, " + b7*", IV_IA2, " + b8*", IV_IA3,
            " + b9*", IV23, CV)
        m.null <- ifelse(is.cubic, paste0(DV, "~ 1 + 0*", IV1,
            " + 0*", IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*",
            IV22, " + 0*", IV13, " + 0*", IV_IA2, " + 0*", IV_IA3,
            " + 0*", IV23, CV), paste0(DV, "~ 1 + 0*", IV1, " + 0*",
            IV2, " + 0*", IV12, " + 0*", IV_IA, " + 0*", IV22,
            CV))
        s.NULL <- lavaan::sem(m.null, data = df[df$out == FALSE, ], fixed.x = TRUE,
            meanstructure = TRUE, se = se, estimator = estimator,
            missing = missing, ...)
    if(!is.null(sampling.weights)){
    df[, sampling.weights] <- data[, sampling.weights]
    sampling.weights <- sampling.weights
}
    if(is.null(sampling.weights)){
    sampling.weights <- NULL
}

    if(!is.null(cluster)){
    df[, cluster] <- data[, cluster]
    cluster <- cluster
}
    if(is.null(cluster)){
    cluster <- NULL
}

    if(!is.null(control.variables)){
    df[, control.variables] <- data[, control.variables]
}




################################
###### Saturated cubic polynomial (full cubic)
        if ("CUBIC" %in% models | is.cubic) {
            if (verbose == TRUE)
                print("Computing full cubic model (CUBIC) ...")
            m.cubic <- paste(polycubic, "u1 := b1 + b2", "u2 := b3 + b4 + b5",
                "u3 := b6 + b7 + b8 + b9", "v1 := b1 - b2", "v2 := b3 - b4 + b5",
                "v3 := b6 - b7 + b8 - b9", add, sep = "\n")
            s.cubic <- lavaan::sem(m.cubic, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }


###### Saturated quadratic polynomial (full quadratic)
if ("QUADRATIC" %in% models ) {
            if (verbose == TRUE)
                print("Computing full quadratic model (QUADRATIC) ...")
            m.quadratic <- paste(polycubic, "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
            "u1 := b1 + b2", "u2 := b3 + b4 + b5",
                "u3 := b6 + b7 + b8 + b9", "v1 := b1 - b2", "v2 := b3 - b4 + b5",
                "v3 := b6 - b7 + b8 - b9", add, sep = "\n")
            s.quadratic <- lavaan::sem(m.quadratic, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }


################################
###### USER SPECIFIC MODEL
          if ("USER" %in% models) {
             list_s.USER <- user_model
               for(i in 1:length(user_model) ){
            m.USER <- paste(polycubic,user_model[[i]],
              	"u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",            
            add, sep = "\n")
            s.USER <- lavaan::sem(m.USER, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
               if (verbose == TRUE)
                print(paste("Computing user specific model (",names(user_model)[i], ")...",sep=""))
            list_s.USER[[ names(user_model)[i] ]] <-    s.USER
                }
              s.USER <- list_s.USER
        }



######################################################################
############### Families of polynomials


################################
###### STEP1 MODELS
    models_STEP1 <- c("FM1_ONLYX","FM2_ONLYY","FM3_ADDITIVE","FM4_INTER","FM5_QUADX","FM6_QUADY","FM7_CONG","FM8_INCONG","FM9_CURVCONGX","FM10_CURVCONGY","FM11_CURVINCONGX","FM12_CURVINCONGY","FM13_QUADXQUADY","FM14_ROTCONG","FM15_ROTINCONG","FM16_CUBICX","FM17_CUBICY","FM18_LEVDEPQUADX","FM19_LEVDEPQUADY","FM20_ASYMCONG","FM21_ASYMINCONG","FM22_LEVDEPCONG","FM23_LEVDEPINCONG","FM24_PARALLELASYM","FM25_NONPARALLELASYM","FM26_PARALLELASYMWEAK","FM27_PARALLELASYMSTRONG","FM28_NONPARALLELASYMWEAK","FM29_NONPARALLELASYMSTRONG","FM30_ASYMCONGROTY","FM31_ASYMCONGROTX","FM32_ASYMINCONGROTY","FM33_ASYMINCONGROTX","FM34_LEVDEPCONGROTY","FM35_LEVDEPCONGROTX","FM36_LEVDEPINCONGROTY","FM37_LEVDEPINCONGROTX")
    if ("STEP1" %in% models) {
        models <- models_STEP1
    }


################################
###### FIRST-ORDER POLYNOMIALS
        # FM1: Single effect X
          if ("FM1_ONLYX" %in% models) {
            if (verbose == TRUE)
                print("Computing only main effect of X model (FM1_ONLYX) ...")
            m.ONLYX <- paste(polycubic, "b2 == 0","b3 == 0","b4 == 0",
                "b5 == 0", "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
                add, sep = "\n")
            s.ONLYX <- lavaan::sem(m.ONLYX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM2: Single effect Y
          if ("FM2_ONLYY" %in% models) {
            if (verbose == TRUE)
                print("Computing only main effect of Y model (FM2_ONLYY) ...")
            m.ONLYY <- paste(polycubic, "b1 == 0","b3 == 0","b4 == 0",
                "b5 == 0", "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ONLYY <- lavaan::sem(m.ONLYY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM3: Additive effects
          if ("FM3_ADDITIVE" %in% models) {
            if (verbose == TRUE)
                print("Computing additive effects model (FM3_ADDITIVE) ...")
            m.ADDITIVE <- paste(polycubic, "b3 == 0","b4 == 0",
                "b5 == 0", "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ADDITIVE <- lavaan::sem(m.ADDITIVE, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

################################
###### SECOND-ORDER POLYNOMIALS
        # FM4: Interaction effect
          if ("FM4_INTER" %in% models) {
            if (verbose == TRUE)
                print("Computing interaction effect model (FM4_INTER) ...")
            m.INTER <- paste(polycubic,"b3 == 0",
                "b5 == 0","b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.INTER <- lavaan::sem(m.INTER, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM5: Quadratic effect X
          if ("FM5_QUADX" %in% models) {
            if (verbose == TRUE)
                print("Computing quadratic effect X model (FM5_QUADX) ...")
            m.QUADX <- paste(polycubic, "b4 == 0","b5 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.QUADX <- lavaan::sem(m.QUADX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM6: Quadratic effect Y
          if ("FM6_QUADY" %in% models) {
            if (verbose == TRUE)
                print("Computing quadratic effect Y model (FM6_QUADY) ...")
            m.QUADY <- paste(polycubic, "b3 == 0","b4 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.QUADY <- lavaan::sem(m.QUADY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

       # FM7: Congruence effect
          if ("FM7_CONG" %in% models) {
            if (verbose == TRUE)
                print("Computing congruence effect model (FM7_CONG) ...")
            m.CONG <- paste(polycubic, "b3 == b4/-2","b3 == b5",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CONG <- lavaan::sem(m.CONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
}

        # FM8: Incongruence effect
          if ("FM8_INCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing incongruence effect model (FM8_INCONG) ...")
            m.INCONG <- paste(polycubic, "b3 == b4/2","b3 == b5",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.INCONG <- lavaan::sem(m.INCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM9: Congruence effect curved by X
          if ("FM9_CURVCONGX" %in% models) {
            if (verbose == TRUE)
                print("Computing congruence effect curved by X model (FM9_CURVCONGX) ...")
            m.CURVCONGX <- paste(polycubic, "b3 == -b4/2","b5 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CURVCONGX <- lavaan::sem(m.CURVCONGX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM10: Congruence effect curved by Y
          if ("FM10_CURVCONGY" %in% models) {
            if (verbose == TRUE)
                print("Computing congruence effect curved by Y model (FM10_CURVCONGY) ...")
            m.CURVCONGY <- paste(polycubic, "b5 == -b4/2","b3 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CURVCONGY <- lavaan::sem(m.CURVCONGY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM11: Incongruence effect curved by X
          if ("FM11_CURVINCONGX" %in% models) {
            if (verbose == TRUE)
                print("Computing incongruence effect curved by X model (FM11_CURVINCONGX) ...")
            m.CURVINCONGX <- paste(polycubic, "b3 == b4/2","b5 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CURVINCONGX <- lavaan::sem(m.CURVINCONGX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM12: Incongruence effect curved by Y
          if ("FM12_CURVINCONGY" %in% models) {
            if (verbose == TRUE)
                print("Computing incongruence effect curved by Y model (FM12_CURVINCONGY) ...")
            m.CURVINCONGY <- paste(polycubic, "b5 == b4/2","b3 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CURVINCONGY <- lavaan::sem(m.CURVINCONGY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM13: Quadratic effects of X and Y
          if ("FM13_QUADXQUADY" %in% models) {
            if (verbose == TRUE)
                print("Computing quadratic effects of X and Y model (FM13_QUADXQUADY) ...")
            m.QUADXQUADY <- paste(polycubic, "b4 == 0",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.QUADXQUADY <- lavaan::sem(m.QUADXQUADY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM14: Rotated congruence effect
          if ("FM14_ROTCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing rotated congruence effect model (FM14_ROTCONG) ...")
            m.ROTCONG <- paste(polycubic, "b4^2 == 4*b3*b5",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ROTCONG <- lavaan::sem(m.ROTCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM15: Rotated incongruence effect
          if ("FM15_ROTINCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing rotated incongruence effect model (FM15_ROTINCONG) ...")
            m.ROTINCONG <- paste(polycubic, "-(b4^2) == 4*b3*b5",
                "b6 == 0", "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ROTINCONG <- lavaan::sem(m.ROTINCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }


################################
###### THIRD-ORDER POLYNOMIALS

        # FM16: Cubic effect of X
          if ("FM16_CUBICX" %in% models) {
            if (verbose == TRUE)
                print("Computing cubic effect of X model (FM16_CUBICX) ...")
            m.CUBICX <- paste(polycubic,
                "b7 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CUBICX <- lavaan::sem(m.CUBICX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM17: Cubic Y effect
          if ("FM17_CUBICY" %in% models) {
            if (verbose == TRUE)
                print("Computing cubic Y effect model (FM17_CUBICY) ...")
            m.CUBICY <- paste(polycubic,
                "b6 == 0", "b7 == 0", "b8 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.CUBICY <- lavaan::sem(m.CUBICY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM18: Level-dependent quadratic effect of X
          if ("FM18_LEVDEPQUADX" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent quadratic effect of X model (FM18_LEVDEPQUADX) ...")
            m.LEVDEPQUADX <- paste(polycubic,
                "b6 == 0", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPQUADX <- lavaan::sem(m.LEVDEPQUADX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM19: Level-dependent quadratic effect of Y
          if ("FM19_LEVDEPQUADY" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent quadratic effect of Y model (FM19_LEVDEPQUADY) ...")
            m.LEVDEPQUADY <- paste(polycubic,
                "b6 == 0", "b7 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPQUADY <- lavaan::sem(m.LEVDEPQUADY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM20: Asymmetric congruence effect
          if ("FM20_ASYMCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric congruence effect model (FM20_ASYMCONG) ...")
            m.ASYMCONG <- paste(polycubic,
                "b6 == b7/-3", "b6 == b8/3", "b6 == -b9",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMCONG <- lavaan::sem(m.ASYMCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM21: Asymmetric incongruence effect
          if ("FM21_ASYMINCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric incongruence effect model (FM21_ASYMINCONG) ...")
            m.ASYMINCONG <- paste(polycubic,
                "b6 == b7/3", "b6 == b8/3", "b6 == b9",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMINCONG <- lavaan::sem(m.ASYMINCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM22: Level-dependent congruence effect
          if ("FM22_LEVDEPCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent congruence effect model (FM22_LEVDEPCONG) ...")
            m.LEVDEPCONG <- paste(polycubic,
                "b6 == -b7", "b6 == -b8", "b6 == b9",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPCONG <- lavaan::sem(m.LEVDEPCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM23: Level-dependent incongruence effect
          if ("FM23_LEVDEPINCONG" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent incongruence effect model (FM23_LEVDEPINCONG) ...")
            m.LEVDEPINCONG <- paste(polycubic,
                "b6 == -b7", "b6 == b8", "b6 == -b9",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPINCONG <- lavaan::sem(m.LEVDEPINCONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM24: Parallel asymmetric congruence and incongruence effects
          if ("FM24_PARALLELASYM" %in% models) {
            if (verbose == TRUE)
                print("Computing parallel asymmetric effects (FM24_PARALLELASYM) ...")
            m.PARALLELASYM <- paste(polycubic,
                "b6 == b8", "b7 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.PARALLELASYM <- lavaan::sem(m.PARALLELASYM, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM25: Non-parallel asymmetric congruence and incongruence effects
          if ("FM25_NONPARALLELASYM" %in% models) {
            if (verbose == TRUE)
                print("Computing non-parallel asymmetric effects (FM25_NONPARALLELASYM) ...")
            m.NONPARALLELASYM <- paste(polycubic,
                "b9 == b7", "b6 == 0", "b8 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.NONPARALLELASYM <- lavaan::sem(m.NONPARALLELASYM, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM26: Parallel asymmetric weak congruence and strong incongruence effects
          if ("FM26_PARALLELASYMWEAK" %in% models) {
            if (verbose == TRUE)
                print("Computing parallel asymmetric weak congruence... model (FM26_PARALLELASYMWEAK) ...")
            m.PARALLELASYMWEAK <- paste(polycubic,
                "b9 == b8/-3", "b6 == 0", "b7 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.PARALLELASYMWEAK <- lavaan::sem(m.PARALLELASYMWEAK, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM27: Parallel asymmetric strong congruence and weak incongruence effects
          if ("FM27_PARALLELASYMSTRONG" %in% models) {
            if (verbose == TRUE)
                print("Computing parallel asymmetric strong congruence... model (FM27_PARALLELASYMSTRONG) ...")
            m.PARALLELASYMSTRONG <- paste(polycubic,
                "b9 == b8/3", "b6 == 0", "b7 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.PARALLELASYMSTRONG <- lavaan::sem(m.PARALLELASYMSTRONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM28: Non-parallel weak congruence and strong incongruence effects
          if ("FM28_NONPARALLELASYMWEAK" %in% models) {
            if (verbose == TRUE)
                print("Computing non-parallel asymmetric weak congruence... model (FM28_NONPARALLELASYMWEAK) ...")
            m.NONPARALLELASYMWEAK <- paste(polycubic,
                "b6 == b7/-3", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.NONPARALLELASYMWEAK <- lavaan::sem(m.NONPARALLELASYMWEAK, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM29: Non-parallel asymmetric strong congruence and weak incongruence effects
          if ("FM29_NONPARALLELASYMSTRONG" %in% models) {
            if (verbose == TRUE)
                print("Computing non-parallel asymmetric strong congruence...model (FM29_NONPARALLELASYMSTRONG) ...")
            m.NONPARALLELASYMSTRONG <- paste(polycubic,
                "b6 == b7/3", "b8 == 0", "b9 == 0",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.NONPARALLELASYMSTRONG <- lavaan::sem(m.NONPARALLELASYMSTRONG, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM30: Asymmetric congruence effect rotated by Y
          if ("FM30_ASYMCONGROTY" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric congruence effect rotated by Y model (FM30_ASYMCONGROTY) ...")
            m.ASYMCONGROTY <- paste(polycubic,
                "b6 == b7/-3", "b6 == b8/3",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMCONGROTY <- lavaan::sem(m.ASYMCONGROTY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM31: Asymmetric congruence effect rotated by X
          if ("FM31_ASYMCONGROTX" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric congruence effect rotated by X model (FM31_ASYMCONGROTX) ...")
            m.ASYMCONGROTX <- paste(polycubic,
                "b9 == b7/3", "b9 == -b8/3",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMCONGROTX <- lavaan::sem(m.ASYMCONGROTX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM32: Asymmetric incongruence effect rotated by Y
          if ("FM32_ASYMINCONGROTY" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric incongruence effect rotated by Y model (FM32_ASYMINCONGROTY) ...")
            m.ASYMINCONGROTY <- paste(polycubic,
                "b6 == b7/3", "b6 == b8/3",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMINCONGROTY <- lavaan::sem(m.ASYMINCONGROTY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM33: Asymmetric incongruence effect rotated by X
          if ("FM33_ASYMINCONGROTX" %in% models) {
            if (verbose == TRUE)
                print("Computing asymmetric incongruence effect rotated by X model (FM33_ASYMINCONGROTX) ...")
            m.ASYMINCONGROTX <- paste(polycubic,
                "b9 == b7/3", "b9 == b8/3",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.ASYMINCONGROTX <- lavaan::sem(m.ASYMINCONGROTX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM34: Level-dependent congruence effect rotated by Y
          if ("FM34_LEVDEPCONGROTY" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent congruence effect rotated by Y model (FM34_LEVDEPCONGROTY) ...")
            m.LEVDEPCONGROTY <- paste(polycubic,
                "b6 == -b7", "b6 == -b8",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPCONGROTY <- lavaan::sem(m.LEVDEPCONGROTY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM35: Level-dependent congruence effect rotated by X
          if ("FM35_LEVDEPCONGROTX" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent congruence effect rotated by X model (FM35_LEVDEPCONGROTX) ...")
            m.LEVDEPCONGROTX <- paste(polycubic,
                "b9 == -b7", "b9 == -b8",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPCONGROTX <- lavaan::sem(m.LEVDEPCONGROTX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM36: Level-dependent incongruence effect rotated by Y
          if ("FM36_LEVDEPINCONGROTY" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent incongruence effect rotated by Y model (FM36_LEVDEPINCONGROTY) ...")
            m.LEVDEPINCONGROTY <- paste(polycubic,
                "b6 == -b7", "b6 == b8",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPINCONGROTY <- lavaan::sem(m.LEVDEPINCONGROTY, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

        # FM37: Level-dependent incongruence effect rotated by X
          if ("FM37_LEVDEPINCONGROTX" %in% models) {
            if (verbose == TRUE)
                print("Computing level-dependent incongruence effect rotated by X model (FM37_LEVDEPINCONGROTX) ...")
            m.LEVDEPINCONGROTX <- paste(polycubic,
                "b9 == -b7", "b9 == b8",
                "u1 := b1 + b2", "u2 := b3 + b4 + b5", "u3 := b6 + b7 + b8 + b9",
                "v1 := b1 - b2", "v2 := b3 - b4 + b5", "v3 := b6 - b7 + b8 - b9",
              add, sep = "\n")
            s.LEVDEPINCONGROTX <- lavaan::sem(m.LEVDEPINCONGROTX, data = df[df$out == FALSE,
                ], fixed.x = TRUE, meanstructure = TRUE, se = se,
                estimator = estimator, missing = missing,sampling.weights=sampling.weights,cluster=cluster, ...)
        }

            })

    modellist <- c(list(CUBIC=s.cubic,
   FM1_ONLYX= s.ONLYX, FM2_ONLYY= s.ONLYY, FM3_ADDITIVE= s.ADDITIVE, FM4_INTER= s.INTER, FM5_QUADX= s.QUADX, FM6_QUADY= s.QUADY,FM7_CONG= s.CONG, FM8_INCONG= s.INCONG, FM9_CURVCONGX= s.CURVCONGX, FM10_CURVCONGY= s.CURVCONGY, FM11_CURVINCONGX= s.CURVINCONGX, FM12_CURVINCONGY= s.CURVINCONGY, FM13_QUADXQUADY= s.QUADXQUADY, FM14_ROTCONG= s.ROTCONG, FM15_ROTINCONG= s.ROTINCONG, FM16_CUBICX= s.CUBICX, FM17_CUBICY= s.CUBICY, FM18_LEVDEPQUADX= s.LEVDEPQUADX, FM19_LEVDEPQUADY= s.LEVDEPQUADY, FM20_ASYMCONG = s.ASYMCONG, FM21_ASYMINCONG= s.ASYMINCONG, FM22_LEVDEPCONG= s.LEVDEPCONG, FM23_LEVDEPINCONG= s.LEVDEPINCONG, FM24_PARALLELASYM= s.PARALLELASYM, FM25_NONPARALLELASYM= s.NONPARALLELASYM, FM26_PARALLELASYMWEAK= s.PARALLELASYMWEAK, FM27_PARALLELASYMSTRONG= s.PARALLELASYMSTRONG, FM28_NONPARALLELASYMWEAK= s.NONPARALLELASYMWEAK, FM29_NONPARALLELASYMSTRONG= s.NONPARALLELASYMSTRONG, FM30_ASYMCONGROTY = s.ASYMCONGROTY, FM31_ASYMCONGROTX = s.ASYMCONGROTX, FM32_ASYMINCONGROTY= s.ASYMINCONGROTY, FM33_ASYMINCONGROTX= s.ASYMINCONGROTX, FM34_LEVDEPCONGROTY= s.LEVDEPCONGROTY, FM35_LEVDEPCONGROTX= s.LEVDEPCONGROTX, FM36_LEVDEPINCONGROTY= s.LEVDEPINCONGROTY, FM37_LEVDEPINCONGROTX= s.LEVDEPINCONGROTX),
   s.USER,QUADRATIC= s.quadratic  )
    res <- list(models = modellist, formula = formula, data = df,
        out.rm = out.rm, outliers = which(df$out == TRUE), DV = DV,
        IV1 = IV1, IV2 = IV2, IV12 = IV12, IV22 = IV22, IV_IA = IV_IA,
        W_IV1 = W_IV1, W_IV2 = W_IV2, IV13 = IV13, IV_IA2 = IV_IA2,
        IV_IA3 = IV_IA3, IV23 = IV23, control.variables = control.variables,
        is.cv = is.cv, r.squared = summary(lm.full)$r.squared,
        is.cubic = is.cubic)
    attr(res, "class") <- "RSA"
    return(res)
}


