#' @title Compare a list of polynomial models against the data
#' @description Compares any number of predefined or user-specific polynomial models and extracts their fit indices, thereby establishing best-fitting solutions.
#' @details This function compares models based on information-theoretic criteria and statistical tests. The cubic saturated polynomial provides a benchmark reference for fit, against which predefined polynomial families (37 to date) or user-specific variants of these families are compared for absolute fit (likelihood ratio test), parsimony (wAIC), explained variance (adjusted R2), and ordinary SEM criteria (e.g., CFI, TLI, RMSEA, SRMR).
#' @param RSA_object x an object of class "RSA_object" generated by RSAmodel()
#' @param order Single or vector of fit indices used to determine best-fitting polynomial families. The output matrix is ordered based on this fit index
#' @param robust Should robust fit indices should be extracted? (default= TRUE)
#' @return A table containing fit indices for each model 
#' 
#' @examples
#' #####ESTIMATE RSA OBJECT
#' RSA_step1 <-  RSAmodel(engagement ~ needs*supplies,
#' data= sim_NSfit, model= c("CUBIC","FM8_INCONG","FM9_INCONG","FM20_ASYMCONG",
#' "FM21_ASYMCONG","FM26_PARALLELASYMWEAK"))
#' ##### COMPARE POLYNOMIAL FAMILIES FROM THE RSA OBJECT
#' RSA_step1_fit <- best.rsa(RSA_step1,order=c("wAIC"))
#' names(RSA_step1$models)
#' #Inspect best-fitting family model
#' summary(RSA_step1$models$FM26_PARALLELASYMWEAK)
#' @export


best.rsa <- function(RSA_object, order=c("wAIC","R2adj"),robust=TRUE){


####Lists to compare
list_models_all <- list()
list_models_fitted <- list()

for(i in 1:length(RSA_object[[1]]) ){
	if(RSA_object[[1]][i]!="NULL" )
list_models_all[[names(RSA_object[[1]][i])]]	<- RSA_object[[1]][i]
	else(NA)
}

for(j in 1:length(list_models_all) ){
		if( is.na(lavaan::lavTestLRT(list_models_all[[j]][[1]])[2,"Chisq"])  )
		NA
		
		else if( !is.na(lavaan::lavTestLRT(list_models_all[[j]][[1]])[2,"Chisq"])  )
list_models_fitted[[names(list_models_all[[j]])]]	<- list_models_all[[j]]
	else(NA)
}


#######Complementary fit indices

#Robust indicators
fit_names_plain	<- c("R2","R2.adj","aic_w","bic_w","aic","bic","df","chisq","pvalue","cfi","tli","rmsea","rmsea.pvalue","srmr")
fit_names_robust	<- c("R2","R2.adj","aic_w","bic_w","aic","bic","df","chisq.scaled","pvalue.scaled","cfi.robust","tli.robust","rmsea.robust","rmsea.pvalue","srmr")	

if(robust==F){ fit_names	<- fit_names_plain }
if(robust==T){ fit_names	<- fit_names_robust	}
names_models <- names(list_models_fitted)
matrix_fitind <- matrix(nrow=length(names_models),ncol=length(fit_names),dimnames=list(names_models,fit_names))
list_sic <- list()

#Extract fit indices
for(m in 1:length(list_models_fitted)){
fitRSA	<- list_models_fitted[[m]][[1]]
matrix_fitind[m,-c(1:4)] <-  lavaan::fitmeasures(fitRSA,fit_names[-c(1:4)])

###R2 and R2adjusted
free.max <- .getFreeParameters(fitRSA)
R2 <- lavaan::inspect(fitRSA, "r2")
matrix_fitind[m,"R2"] <- R2
N <- lavaan::nobs(fitRSA)
k <- free.max - matrix_fitind[m,"df"]
matrix_fitind[m,"R2.adj"] <- 1 - ((1 - R2)) * ((N - 1)/(N - k - 1))
}
colnames(matrix_fitind) <- fit_names_plain


#Akaike and bayesian weights
delta_i_aic <- matrix_fitind[, "aic"] - min(matrix_fitind[, "aic"],na.rm=T)
matrix_fitind[,"aic_w"] <- exp(-1/2* delta_i_aic)/sum(exp(-1/2* delta_i_aic))
delta_i_bic <- matrix_fitind[, "bic"] - min(matrix_fitind[, "bic"],na.rm=T)
matrix_fitind[,"bic_w"] <- exp(-1/2* delta_i_bic)/sum(exp(-1/2* delta_i_bic))
matrix_fitind <- round(data.frame(matrix_fitind),3)
head(matrix_fitind)


#### Merge fit indices
merge_fit <- cbind(matrix_fitind[,c("chisq","df","pvalue")], matrix_fitind[,c("aic_w","bic_w","aic","bic")], R2= matrix_fitind[,c("R2")],R2adj= matrix_fitind[,c("R2.adj")], matrix_fitind[,c("cfi","tli","rmsea","rmsea.pvalue","srmr")])
merge_fit <- round(data.frame(merge_fit),3)
merge_fit[,c("chisq","aic","bic")] <- round(merge_fit[,c("chisq","aic","bic")],1)
head(merge_fit)

#### Fit indices out
# Column names
head(merge_fit)
colnames(merge_fit) <- c("LRT_chi2","LRT_df","LRT_pvalue","wAIC","wBIC","AIC","BIC","R2","R2adj","CFI","TLI","RMSEA","RMSEA_pvalue","SRMR")

# Ordered by "order"
merge_fit_order <- merge_fit[-1,]
order_all <- merge_fit_order[, order[1]]
if(length(order)>1){
order_all <- list()	
for(o in 1:length(order)){
order_all[[ order[o] ]]  <- merge_fit_order[, order[o]]
}
}

if(length(order)==1){
rsa_fit_order <- merge_fit_order[order(order_all,decreasing=T,na.last=F),]
}

if(length(order)==2){
rsa_fit_order <- merge_fit_order[order(order_all[[1]],order_all[[2]],decreasing=T,na.last=F),]
}

if(length(order)==3){
rsa_fit_order <- merge_fit_order[order(order_all[[1]],order_all[[2]],order_all[[3]],decreasing=T,na.last=F),]
}
rsa_fit_out <- data.frame(rbind(merge_fit[1,], rsa_fit_order))
head(rsa_fit_out)



##### OUT all
compare_out <- rsa_fit_out[,c("LRT_chi2","LRT_df","LRT_pvalue","AIC","wAIC","R2adj","CFI","TLI","RMSEA","SRMR","BIC","wBIC","R2")]  
compare_out

}
