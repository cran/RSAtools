#' @title Plot extrema in the response surface along the lines of congruence and incongruence
#'
#' @description
#' Plots the response surface of an RSA object along the lines of congruence (LOC) and incongruence (LOIC), while indicating the location of reversal or acceleration points in the surface (when they exist).
#'
#' @details
#' The lines of congruence (LOC) and incongruence (LOIC) are fundamental to RSA applications. This function allows plotting the response along the LOC and LOIC, thus offering a more precise visualization of the response surface (as a complement to 3d or 2d projections provided by \code{RSAplot}). The location of reversal or acceleration points can be provided with x-y-z coordinales or without them (e_label='none').
#' @param RSA_object An RSA object
#' @param model The model to be probed for extrema (reversal or acceleration points)
#' @param acceleration Rates of accelerations along the LOC and LOIC to be inspected (0< |rate| < 1). Passed on internally to \code{ident.ext}
#' @param n_sample Number of random draws to consider to find extrema. This option is used for large samples to increase speed in preliminary analyses, but it is not recommended for published results). Passed on internally to \code{ident.ext}  
#' @param names_xLOC Label of x axis on the LOC plot.
#' @param names_xLOIC Label of x axis on the LOIC plot.
#' @param names_z Label of z axis (outcome).
#' @param xlim Limits of the x axis
#' @param zlim Limits of the z axis
#' @param e_label If "none", no coordinates are projected. Defaults to NULL.
#' @param text_size Text size for titles and labels.
#' @param elabel_size Text size for extrema coordinates.
#' @param e_size Point size of extrema.
#' @return A list of plot of the lines of congruence (LOC) and incongruence (LOIC)
#' @seealso \code{\link{RSAmodel}}
#'
#' @examples
#' ######PLOT EXTREMA OVER LOC AND LOIC
#' EXTsim <- plotting.ext(RSA_step1,model="FM26_PARALLELASYMWEAK",
#' xlim=c(-3,3),zlim=c(-3,3),acceleration=c(0,-0.3),
#' text_size=0.7,elabel_size=3,e_size=2)
#' ggpubr::ggarrange(EXTsim[["LOC"]], EXTsim[["LOIC"]],
#' labels = c("A. Response over LOC", "B. Response over LOIC"),
#' nrow=1,ncol=2,font.label = list(size = 11))
#' @export


plotting.ext <- function(RSA_object, model,acceleration=c(0,0),n_sample= 100,names_xLOC=NULL, names_xLOIC=NULL, names_z=NULL,xlim=NULL,zlim=NULL, e_label=NULL, text_size=1,elabel_size=5, e_size=3){

###var_x, var_z
			if (is.null(names_xLOC)) {
				name_x <- RSA_object$IV1
				name_y <- RSA_object$IV2
				names_xLOC <- paste(name_x," = ", name_y,sep="")
				}
			if (is.null(names_z)) {
				names_z <- RSA_object$DV }

			if (is.null(names_xLOIC)) {
				names_xLOIC <- paste(name_x," = -", name_y,sep="")
				 }

			
###xlim, zlim
data = RSA_object$data
			if (is.null(xlim)) {
				# xlim <- c(-3, 3)
				xlim <- c(min(data[,c(2,3)],na.rm=T),max(data[,c(2,3)],na.rm=T))
				# xlim[1] <- xlim[1]*ifelse(xlim[1]>0, 1.2, 0.8)
				# xlim[2] <- xlim[2]*ifelse(xlim[2]>0, 0.8, 1.2)
				}
			if (is.null(zlim)) {
				# zlim <- c(-3, 3)
				zlim <- c(min(data[,c(1)],na.rm=T),max(data[,c(1)],na.rm=T))
				# zlim[1] <- zlim[1]*ifelse(zlim[1]<0, 1.1, 0.9)
				# zlim[2] <- zlim[2]*ifelse(zlim[2]<0, 0.9, 1.1)
				}			

#Data range and model
data = RSA_object$data
names(data)
# x_min <-  min(data[,2],na.rm=T)
# x_max <-  max(data[,c(2)],na.rm=T)
xlim_raw <- xlim
var_x <- seq(xlim_raw[1], xlim_raw[2], psych::describe(seq(xlim_raw[1], xlim_raw[2]))$se/100)
C <- round(coef(RSA_object$models[[model]]),3)


###Get intercept value
####Adjust intercept for control variables (if any)
b0 = C[paste(RSA_object$DV,"~1",sep="")]
 	if(RSA_object$is.cv){
	  # vector of means of the control variables
	  cvmeans <- colMeans(RSA_object$data[, RSA_object$control.variables], na.rm=T)
	  # coefficients of the control variables
	  cvcoefs <- C[paste0(RSA_object$DV, "~", RSA_object$control.variables)]
	  # new b0 = b0 + cv1*mean(controlvariable1) + cv2*mean(controlvariable2) + ...
	  b0 <- b0 + sum(cvmeans*cvcoefs)	  
	  }


 #Df LOC and LOIC
u1 <- ( C["b1"]+ C["b2"] )
u2 <- ( C["b3"]+C["b4"]+C["b5"] )
u3 <- (C["b6"]+C["b7"]+C["b8"]+C["b9"])
z_LOC <- b0+u1*var_x + u2* var_x^2 + u3* var_x^3 
df_LOC <- data.frame(var_x, z_LOC)
v1 <- ( C["b1"]- C["b2"])
v2 <- ( C["b3"]-C["b4"]+C["b5"] )
v3 <- (C["b6"]-C["b7"]+C["b8"]-C["b9"])
z_LOIC <- b0+v1*var_x + v2* var_x^2 + v3* var_x^3 
df_LOIC <- data.frame(var_x, z_LOIC)

 #Compute extrema
extrema <- ident.ext(RSA_object,alpha=0.05,model=model, n_sample= n_sample,acceleration= acceleration)
extrema

 #Acceleration points LOC				
acceleration_loc <- extrema[c("a1_LOC","a2_LOC"),c("X_value","Y_value","Z_value")]
acceleration_loc <- acceleration_loc[which(acceleration_loc$X_value!="no ext."),]
acceleration_loc$Z_value <- as.numeric(acceleration_loc$Z_value)
acceleration_loc$X_value <- as.numeric(acceleration_loc$X_value)
acceleration_loc$Y_value <- as.numeric(acceleration_loc$Y_value)

if(length(acceleration_loc$X_value)!=0){
label_r1_LOC <- ifelse(!is.na(acceleration_loc["a1_LOC","X_value"]),paste("a1_LOC: ","\n","X = ", acceleration_loc[c("a1_LOC"),"X_value"],"; Z = ", acceleration_loc[c("a1_LOC"),"Z_value"],sep=""),"r1_LOC/a1_LOC: none")
label_r2_LOC <- ifelse(!is.na(acceleration_loc["a2_LOC","X_value"]),paste("a2_LOC: ","\n","X = ", acceleration_loc[c("a2_LOC"),"X_value"],"; Z = ", acceleration_loc[c("a2_LOC"),"Z_value"],sep=""),"r2_LOC/a2_LOC: none")
labels_loc <- paste(label_r1_LOC, label_r2_LOC,sep="\n")
}


 #Reversal points LOC				
reversal_loc <- extrema[c("r1_LOC","r2_LOC"),c("X_value","Y_value","Z_value")]
reversal_loc <- reversal_loc[which(reversal_loc$X_value!="no ext."),]
reversal_loc$Z_value <- as.numeric(reversal_loc$Z_value)
reversal_loc$X_value <- as.numeric(reversal_loc$X_value)
reversal_loc$Y_value <- as.numeric(reversal_loc$Y_value)

if(length(acceleration_loc$X_value)==0  ){
label_r1_LOC <- ifelse(!is.na(reversal_loc["r1_LOC","X_value"]),paste("r1_LOC: ","\n","X = ", reversal_loc[c("r1_LOC"),"X_value"],"; Z = ", reversal_loc[c("r1_LOC"),"Z_value"],sep=""),"r1_LOC/a1_LOC: none")
label_r2_LOC <- ifelse(!is.na(reversal_loc["r2_LOC","X_value"]),paste("r2_LOC: ","\n","X = ", reversal_loc[c("r2_LOC"),"X_value"],"; Z = ", reversal_loc[c("r2_LOC"),"Z_value"],sep=""),"r2_LOC/a2_LOC: none")
labels_loc <- label_r1_LOC
}


 #Acceleration points LOIC				
acceleration_loic <- extrema[c("a1_LOIC","a2_LOIC"),c("X_value","Y_value","Z_value")]
acceleration_loic <- acceleration_loic[which(acceleration_loic$X_value!="no ext."),]
acceleration_loic$Z_value <- as.numeric(acceleration_loic$Z_value)
acceleration_loic$X_value <- as.numeric(acceleration_loic$X_value)
acceleration_loic$Y_value <- as.numeric(acceleration_loic$Y_value)

if(length(acceleration_loic$X_value)!=0){
label_r1_LOIC <- ifelse(!is.na(acceleration_loic["a1_LOIC","X_value"]),paste("a1_LOIC: ","\n","X = ", acceleration_loic[c("a1_LOIC"),"X_value"],"; Z = ", acceleration_loic[c("a1_LOIC"),"Z_value"],sep=""),"r1_LOIC/a1_LOIC: none")
label_r2_LOIC <- ifelse(!is.na(acceleration_loic["a2_LOIC","X_value"]),paste("a2_LOIC: ","\n","X = ", acceleration_loic[c("a2_LOIC"),"X_value"],"; Z = ", acceleration_loic[c("a2_LOIC"),"Z_value"],sep=""),"r2_LOIC/a2_LOIC: none")
labels_LOIC <- paste(label_r1_LOIC, label_r2_LOIC,sep="\n")
}


 #Reversal points LOIC				
reversal_loic <- extrema[c("r1_LOIC","r2_LOIC"),c("X_value","Y_value","Z_value")]
reversal_loic <- reversal_loic[which(reversal_loic$X_value!="no ext."),]
reversal_loic$Z_value <- as.numeric(reversal_loic$Z_value)
reversal_loic$X_value <- as.numeric(reversal_loic$X_value)
reversal_loic$Y_value <- as.numeric(reversal_loic$Y_value)

if(length(acceleration_loic$X_value)==0){
label_r1_LOIC <- ifelse(!is.na(reversal_loic["r1_LOIC","X_value"]),paste("r1_LOIC: ","\n","X = ", reversal_loic[c("r1_LOIC"),"X_value"],"; Z = ", reversal_loic[c("r1_LOIC"),"Z_value"],sep=""),"r1_LOIC/a1_LOIC: none")
label_r2_LOIC <- ifelse(!is.na(reversal_loic["r2_LOIC","X_value"]),paste("r2_LOIC: ","\n","X = ", reversal_loic[c("r2_LOIC"),"X_value"],"; Z = ", reversal_loic[c("r2_LOIC"),"Z_value"],sep=""),"r2_LOIC/a2_LOIC: none")
labels_loic <- label_r1_LOIC
}



####CANCEL PROJECTION OF EXTREMA LABELS
if(!is.null(e_label) ){
label_r1_LOC <- NA
label_r2_LOC <- NA
label_r1_LOIC <- NA
label_r2_LOIC <- NA
}


####EXTREMA COORDINATES
xLOC_e1 <- ifelse(is.na(reversal_loc[1,"X_value"]),acceleration_loc[1,"X_value"],reversal_loc[1,"X_value"])
zLOC_e1 <- ifelse(is.na(reversal_loc[1,"Z_value"]), acceleration_loc[1,"Z_value"],reversal_loc[1,"Z_value"])
xLOC_e2 <- ifelse(is.na(reversal_loc[2,"X_value"]),acceleration_loc[2,"X_value"],reversal_loc[2,"X_value"])
zLOC_e2 <- ifelse(is.na(reversal_loc[2,"Z_value"]), acceleration_loc[2,"Z_value"],reversal_loc[2,"Z_value"])
xLOIC_e1 <- ifelse(is.na(reversal_loic[1,"X_value"]),acceleration_loic[1,"X_value"],reversal_loic[1,"X_value"])
zLOIC_e1 <- ifelse(is.na(reversal_loic[1,"Z_value"]), acceleration_loic[1,"Z_value"],reversal_loic[1,"Z_value"])
xLOIC_e2 <- ifelse(is.na(reversal_loic[2,"X_value"]),acceleration_loic[2,"X_value"],reversal_loic[2,"X_value"])
zLOIC_e2 <- ifelse(is.na(reversal_loic[2,"Z_value"]), acceleration_loic[2,"Z_value"],reversal_loic[2,"Z_value"])


#####PLOT SOLUTIONS
code_color <- c("black","black","red","green","green")
code_line <- c("solid","twodash","twodash")

####PLOT Z_LOC
df_LOC$groupEXT <- ifelse(df_LOC$var_x <xLOC_e1 & !is.na(xLOC_e2),"2_ext1",ifelse(df_LOC$var_x >xLOC_e2 & !is.na(xLOC_e2),"3_ext1",ifelse(is.na(xLOC_e2),"Quadratic",ifelse(is.na(xLOC_e2) & is.na(xLOC_e1),"None","1_Inflection"))))
groupEXT <- df_LOC$groupEXT

if( levels(factor(df_LOC$groupEXT))[1]=="None"  ){
plot_zLOC <- ggplot(df_LOC,aes(x= var_x,y= z_LOC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(mean(var_x), mean(z_LOC),label = paste(label_r1_LOC,"\n", label_r2_LOC,sep=""), vjust = 0),size=elabel_size,colour="green",inherit.aes=F, fontface = "plain")+geom_point(aes(mean(var_x), mean(z_LOC)), size=e_size, shape="circle", color="green")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}


if( levels(factor(df_LOC$groupEXT))[1]=="1_Inflection"  ){
plot_zLOC <- ggplot(df_LOC,aes(x= var_x,y= z_LOC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
# geom_text(aes(mean(var_x), mean(z_LOC),label = paste(label_r1_LOC,"\n", label_r2_LOC,sep=""), vjust = 0),size=elabel_size,colour="green",inherit.aes=F, fontface = "plain")+geom_point(aes(mean(var_x), mean(z_LOC)), size=e_size, shape="circle", color="green")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}



if( levels(factor(df_LOC$groupEXT))[1]=="Quadratic" | length(table(df_LOC$groupEXT))==2 ){
plot_zLOC <- ggplot(df_LOC,aes(x= var_x,y= z_LOC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(xLOC_e1, zLOC_e1,label = label_r1_LOC, vjust = -u2/abs(u2)* -0.5),size=elabel_size,colour="green",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOC_e1, zLOC_e1), size=e_size, shape="circle", color="green")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}


code_line <- c("twodash","solid","solid")
if( length(table(df_LOC$groupEXT))>2 ){
plot_zLOC <- ggplot(df_LOC,aes(x= var_x,y= z_LOC,group=groupEXT,linetype= groupEXT) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(xLOC_e1, zLOC_e1,label = label_r1_LOC, vjust = xLOC_e1/abs(xLOC_e1)*0.8),size=elabel_size,colour="green",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOC_e1, zLOC_e1), size=e_size, shape="circle", color="green")+
geom_text(aes(xLOC_e2, zLOC_e2,label = label_r2_LOC, vjust = xLOC_e2/abs(xLOC_e2)*1.8),size=elabel_size,colour="green",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOC_e2, zLOC_e2), size=e_size, shape="circle", color="green")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}



####PLOT Z_LOIC
code_line <- c("solid","twodash","twodash")
df_LOIC$groupEXT2 <- ifelse(df_LOIC$var_x < xLOIC_e1 & !is.na(xLOIC_e2),"2_ext1",ifelse(df_LOIC$var_x >xLOIC_e2 & !is.na(xLOIC_e2),"3_ext1",ifelse(is.na(xLOIC_e2),"Quadratic",ifelse(is.na(xLOIC_e2) & is.na(xLOIC_e1),"None","1_Inflection"))))
groupEXT2 <- df_LOIC$groupEXT2

if( levels(factor(df_LOIC$groupEXT2))[1]=="None"  ){
plot_zLOIC <- ggplot(df_LOIC,aes(x= var_x,y= z_LOIC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOIC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(mean(var_x), mean(z_LOIC),label = paste(label_r1_LOIC,"\n", label_r2_LOIC,sep=""), vjust = 0),size=elabel_size,colour="red",inherit.aes=F, fontface = "plain")+geom_point(aes(mean(var_x), mean(z_LOIC)), size=e_size, shape="circle", color="red")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}

if(  levels(factor(df_LOIC$groupEXT2))[1]=="1_Inflection" ){
plot_zLOIC <- ggplot(df_LOIC,aes(x= var_x,y= z_LOIC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x=names_xLOIC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
# geom_text(aes(mean(var_x), mean(z_LOIC),label = paste(label_r1_LOIC,"\n", label_r2_LOIC,sep=""), vjust = 0),size=elabel_size,colour="red",inherit.aes=F, fontface = "plain")+geom_point(aes(mean(var_x), mean(z_LOIC)), size=e_size, shape="circle", color="red")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}


if( levels(factor(df_LOIC$groupEXT2))[1]=="Quadratic" | length(table(df_LOIC$groupEXT2))==2){
plot_zLOIC <- ggplot(df_LOIC,aes(x= var_x,y= z_LOIC) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x= names_xLOIC,y= names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(xLOIC_e1, zLOIC_e1,label = label_r1_LOIC, vjust = -v2/abs(v2)* -0.5),size=elabel_size,colour="red",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOIC_e1, zLOIC_e1), size=e_size, shape="circle", color="red")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}

code_line <- c("twodash","solid","solid")
if( length(table(df_LOIC$groupEXT2))>2    ){
plot_zLOIC <- ggplot(df_LOIC,aes(x= var_x,y= z_LOIC,group=groupEXT2,linetype= groupEXT2) )+ggtitle("")+geom_line(na.rm=TRUE,size=1)+
labs(x= names_xLOIC,y=names_z)+
scale_linetype_manual(values= code_line)+scale_color_manual(values=code_color)+
scale_y_continuous(breaks=seq(-10, 10,0.5),limits=zlim)+
scale_x_continuous(breaks=seq(-10, 10,0.5),limits=xlim)+
geom_text(aes(xLOIC_e1, zLOIC_e1,label = label_r1_LOIC, vjust = xLOIC_e1/abs(xLOIC_e1)*0.8),size=elabel_size,colour="red",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOIC_e1, zLOIC_e1), size=e_size, shape="circle", color="red")+
geom_text(aes(xLOIC_e2, zLOIC_e2,label = label_r2_LOIC, vjust = xLOIC_e2/abs(xLOIC_e2)*1.8),size=elabel_size,colour="red",inherit.aes=F, fontface = "plain")+geom_point(aes(xLOIC_e2, zLOIC_e2), size=e_size, shape="circle", color="red")+
theme(legend.position="none",panel.background=element_rect(fill="grey100",colour="grey50"),legend.key.size = unit(2, units="cm"),legend.background=element_rect(colour="black"),panel.grid=element_line(color="grey"),axis.title=element_text(size=15*text_size), axis.text=element_text(size=15*text_size));

}

list_plot <- list(LOC=plot_zLOC, LOIC=plot_zLOIC) 
return(list_plot)

}