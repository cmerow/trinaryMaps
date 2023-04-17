#################################################################
#################################################################
#################################################################
#' Plot trinary map raster
#'
#' This function plots the trinary map raster for a given species. The raster will have three colors: grey for the area with no prediction, blue for the predicted presence area, and red for the predicted absence area. Additionally, if the \code{pres} argument is provided, the function will overlay presence data in grey crosses on the plot. If the \code{shapesToPlot} argument is provided, the function will overlay the specified shapefiles on the plot. If a \code{plotFile} is provided, the plot will be saved to that file. If \code{openFig} is TRUE, the plot will also be opened in the system's default PDF viewer.
#'
#' @param trinaryRaster Trinary map raster created by \code{\link{trinaryMapWorkflow}}.
#' @param plotFile File to save plot. If NULL, plot is not saved.
#' @param pres Presence data to overlay on plot. If NULL, presence data is not included.
#' @param main Title for plot. Default is NULL.
#' @param shapesToPlot List of shapefiles to overlay on plot. Default is NULL.
#' @param openFig Logical indicating whether to open plot in system's default PDF viewer. Default is TRUE.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' # plot trinary map for species 1
#' trinaryMapPlot(trinaryRaster = trinaryMapSDMPW1$trinaryMap,
#' pres = species1Presences,
#' shapesToPlot = list(world.shp),
#' plotFile = "trinary_map_species1.pdf",
#' openFig = TRUE)
#' }

trinaryMapPlot=function(trinaryRaster,
												plotFile=NULL,
												pres=NULL,
												main=NULL,					
												shapesToPlot=NULL,
												openFig=TRUE){

	#  for testing
	#  
	if(!is.null(plotFile))  grDevices::pdf(plotFile,h=7,w=7)

	if(!is.null(pres)){
		if(is.na(projection(pres))) {
			projection(pres)=projection(trinaryRaster)
			message('The projection for your presences was not specified, so assuming its the same as your rasters')
		}
	}
	if(!is.null(plotFile))  grDevices::pdf(plotFile,h=7,w=7)

		raster::image(trinaryRaster==0, xlab='',ylab='',col=c(grey(1,0),'grey80'),xaxt='n', yaxt='n',bty='n',main=main)
		raster::image(trinaryRaster==1, xlab='',ylab='',col=c(grey(1,0),'steelblue'),xaxt='n', yaxt='n',bty='n',add=TRUE)
		raster::image(trinaryRaster>=2, xlab='',ylab='',col=c(grey(1,0),'red1'),xaxt='n', yaxt='n',bty='n',add=TRUE)


		if(!is.null(pres)) graphics::points(pres,col='grey90',pch='+',cex=.4)

		if(!is.null(shapesToPlot)) lapply(shapesToPlot,sp::plot,add=T,lwd=.1,border='grey40')

	if(!is.null(plotFile))  grDevices::dev.off(); if(openFig) system(paste0('open ',plotFile))
	
}


#################################################################
#################################################################
#################################################################
#' Plot trinary ROC curve and save the output as pdf file.
#'
#' @param trinaryPlotThings Object of trinary ROC curve, contains xx, y, xout, and y. as its attributes.
#' @param trinaryDF Data frame with trinary ROC curve information, including youden.thresh.roc.x, hi.thresh.roc.x, lo.thresh.roc.x, and trinary.pauc.
#' @param plotFile Path to save the output as pdf file. If NULL, the plot will not be saved.
#' @param openFig Logical, if TRUE, the plot will be shown in a new window.
#' @return None

#' @export

trinaryROCPlot=function(trinaryPlotThings,
												trinaryDF,
												plotFile=NULL,
												openFig=TRUE
												#addMapPlot=FALSE
												){
	#  for testing
	#  trinaryPlotThings=threshs$trinaryPlotThings; trinaryDF=threshs[[1]]; openFig=T; 												addMapPlot=FALSE

	if(!is.null(plotFile))  grDevices::pdf(plotFile,h=1.6*7,w=7.6)
		par(mfrow=c(3,2),mar=c(5,4,2,2))
		plot(trinaryPlotThings$xx,trinaryPlotThings$y,type='l',lwd=4,ylim=c(0,1.08), main="ROC",xlab='1-specificity', ylab="sensitivity",xlim=c(-.08,1))
		#abline(h=0,lty=2)
		abline(v=trinaryDF$youden.thresh.roc.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=trinaryDF$hi.thresh.roc.x,col='steelblue',lty=3,lwd=2)
		abline(v=trinaryDF$lo.thresh.roc.x,col='red1',lty=3,lwd=2)

		abline(0,1,lty=2)
		#points(trinaryDF$lo.thresh.roc.x,y.lo,pch=21,bg='red1')
		points(trinaryDF$youden.thresh.roc.x,trinaryDF$youden.thresh.roc.y,pch=21, bg='darkgoldenrod2')
		points(trinaryDF$hi.thresh.roc.x,trinaryDF$hi.thresh.roc.y,pch=21, bg='steelblue')
		points(trinaryDF$lo.thresh.roc.x,trinaryDF$lo.thresh.roc.y,pch=21,bg='red1')
		text(.8,.2,paste0('pAUC =\n',round(trinaryDF$trinary.pauc,3)),cex=1.5)
		text(trinaryDF$youden.thresh.roc.x+.02,trinaryDF$youden.thresh.roc.y-.06, paste0('(',round(trinaryDF$youden.thresh.roc.x,2),',', round(trinaryDF$youden.thresh.roc.y,2),')'), adj=c(0,0), cex=1,col='darkgoldenrod2') 
		text(trinaryDF$hi.thresh.roc.x+.02,trinaryDF$hi.thresh.roc.y+.04, paste0('(',round(trinaryDF$hi.thresh.roc.x,2),',',round(trinaryDF$hi.thresh.roc.y,2),')'),adj=c(0,0), cex=1,col='steelblue')
		text(trinaryDF$lo.thresh.roc.x+.02,trinaryDF$lo.thresh.roc.y-.06, paste0('(',round(trinaryDF$lo.thresh.roc.x,2),',',round(trinaryDF$lo.thresh.roc.y,2),')'),adj=c(0,1), cex=1,col='red1')
		
		# ROC' -----------------------------------------
		if(is.null(trinaryPlotThings$y.)) {
				plot(.5,.5,col='white', main="ROC'", xlab='1-specificity',ylab="sensitivity'");
				graphics::text(.5,.5,'Derivatives were not used for this model.\n 0% and 30% training presence quantiles\n were used instead')
		} else {
			if(any(is.na(trinaryPlotThings$y.))){ #this catches NAs 
				keep=which(!is.na(trinaryPlotThings$y.))
				x.tmp=trinaryPlotThings$xout[keep]; y.tmp=trinaryPlotThings$y.[keep]
			} else { x.tmp=trinaryPlotThings$xout; y.tmp=trinaryPlotThings$y. }
				# catch infs; don't need this for second derivs, because those become NA
			keep=which(is.finite(y.tmp))
			y.tmp=y.tmp[keep]; x.tmp=x.tmp[keep]
			if(any(is.infinite(y.tmp))){
				plot(.5,.5,col='white', main="ROC'", xlab='1-specificity',ylab="sensitivity'"); graphics::text(.5,.5,'Derivative contains Inf')
			} else {
					p1=plot(x.tmp,y.tmp,type='l',col='black',log='y', main="ROC'", xlab='1-specificity',ylab="sensitivity''",xlim=c(0,1))
			}
	
			abline(h=0,lty=2)
			abline(h=1,lty=2)
			abline(v=trinaryDF$youden.thresh.roc.x,col='darkgoldenrod2',lty=3,lwd=2)
			abline(v=trinaryDF$hi.thresh.roc.x,col='steelblue',lty=3,lwd=2)
			abline(v=trinaryDF$lo.thresh.roc.x,col='red1',lty=3,lwd=2)
		}
		
		# ROC'' -----------------------------------------
		if(is.null(trinaryPlotThings$y..)) {
				plot(.5,.5,col='white', main="logmod(ROC'')", xlab='1-specificity',ylab="sensitivity'");
				graphics::text(.5,.5,'Derivatives were not used for this model.\n 0% and 30% training presence quantiles\n were used instead')
		} else {
			if(any(is.na(trinaryPlotThings$y..))){
				keep=which(!is.na(trinaryPlotThings$y..))
				x.tmp=trinaryPlotThings$xout[keep]; y.tmp=trinaryPlotThings$y..[keep]
			} else {x.tmp=trinaryPlotThings$xout; y.tmp=trinaryPlotThings$y..}
			if(any(is.infinite(y.tmp)) | all(y.tmp==0)){
				plot(.5,.5,col='white', main="logmod(ROC'')", xlab='1-specificity',ylab="logmod(sensitivity'')"); graphics::text(.5,.5,'Derivative contains Inf')
			} else {
				plot(x.tmp,logmod(y.tmp),type='l',col='black', main="logmod(ROC'')", xlim=c(0,1),xlab='1-specificity',ylab="logmod(sensitivity'')")
			}

			abline(h=0,lty=2) 
			abline(v=trinaryDF$youden.thresh.roc.x,col='darkgoldenrod2',lty=3,lwd=2)
			abline(v=trinaryDF$hi.thresh.roc.x,col='steelblue',lty=3,lwd=2)
			abline(v=trinaryDF$lo.thresh.roc.x,col='red1',lty=3,lwd=2)
		}
		
		# inverse ROC -----------------------------------
		if(is.null(trinaryPlotThings$y1)) {
				plot(.5,.5,col='white', main="Inverse ROC'", xlab='1-sensitivity',ylab='specificity');
				graphics::text(.5,.5,'Derivatives were not used for this model.\n 0% and 30% training presence quantiles\n were used instead')
		} else {
			if(any(is.na(trinaryPlotThings$y1))){
				keep=which(!is.na(trinaryPlotThings$y1))
				x.tmp=trinaryPlotThings$xx1[keep]; y.tmp=trinaryPlotThings$y1[keep]
			} else {x.tmp=trinaryPlotThings$xx1; y.tmp=trinaryPlotThings$y1}
			if(any(is.infinite(y.tmp))){
				plot(.5,.5,col='white', main="Inverse ROC", xlab='1-sensitivity',ylab='specificity'); graphics::text(.5,.5,'Derivative contains Inf')
			} else {
					plot(x.tmp,y.tmp,type='l',lwd=4,ylim=c(0,1.08), main='Inverse ROC',xlab='1-sensitivity',ylab='specificity')
			}

			abline(h=trinaryDF$y.lo.inv,col='red1',lty=3,lwd=2)
			abline(v=trinaryDF$x.lo.inv,col='red1',lty=3,lwd=2)
			abline(0,1,lty=2)
		}
		
		# inverse ROC'' -----------------------------------
		if(is.null(trinaryPlotThings$y1..)) {
				plot(.5,.5,col='white',main="logmod[(Inverse ROC)'']", xlab='1-sensitivity',ylab="logmod(specificity'')");
				graphics::text(.5,.5,'Derivatives were not used for this model.\n 0% and 30% training presence quantiles\n were used instead')
		} else {
			if(!any(is.finite(trinaryPlotThings$y1..))){
				keep=which(!is.na(trinaryPlotThings$y1..))
				x.tmp=trinaryPlotThings$x1out[keep]; y.tmp=trinaryPlotThings$y1..[keep]
			} else {x.tmp=trinaryPlotThings$x1out; y.tmp=trinaryPlotThings$y1..}
			if(any(is.infinite(y.tmp)) | all(y.tmp==0,na.rm=T)){
				plot(.5,.5,col='white', main="logmod[(Inverse ROC)'']", xlab='1-sensitivity',ylab="logmod(specificity'')"); graphics::text(.5,.5,'Derivative contains Inf')
			} else {
				plot(x.tmp,logmod(y.tmp),type='l', main="logmod[(Inverse ROC)'']",xlab='1-sensitivity',ylab="logmod(specificity'')",xlim=c(0,1),ylim=range(logmod(trinaryPlotThings$y1..),na.rm=T))
			}
		
			#abline(h=y.lo) # not sure what this was supposed to be... was defined as 1-x.lo.inv in the modeling....
			abline(h=0,lty=2)
			abline(v=trinaryDF$x.lo.inv,col='red1',lty=3,lwd=2)
		}
		# for third derivative - works but not needed
		# # 		plot(trinaryPlotThings$xout,logmod(y...),type='l',col='black', log='y',xlim=c(0,1),main="ROC'''")
		# # 		abline(h=0,lty=2)
		# # 		abline(v=x.lo,col='darkgoldenrod2',lty=3,lwd=2)
		# # 		abline(v=trinaryDF$hi.thresh.roc.x,col='steelblue',lty=3,lwd=2)
		# # 		abline(v=x.lo,col='red1',lty=3,lwd=2)

		# curvature; code works, idea doesn't
		# # 		plot(trinaryPlotThings$xout,c2,type='l',col='black',xlim=c(0,1), main=expression(paste(kappa,"(ROC)")))
		# # 		abline(v=x.lo,col='red1',lty=3,lwd=2)
		# # 		lines(trinaryPlotThings$xout,c2.1,col='grey30')
		# # 		lines(trinaryPlotThings$xout,c2.2,col='grey60')
		#abline(v=.y..f(rate,skew),lty=2)
		##abline(v=x.root,col='blue',lty=2)
		#lines(ss$x,ss$y,col='blue')
		
# 		if(addMapPlot){
# 			print('not done')
# 			# have to push all these args, which is annoying
# 			# trinaryMapPlot(trinaryRaster,plotFile,pres=pres, expertRasterPath=expertRasterPath,expertShpPath=expertShpPath,shapesToPlot=shapesToPlot,openFig=T)
# 		}

	
	if(!is.null(plotFile))  grDevices::dev.off()
	if(openFig & !is.null(plotFile))	system(paste0('open ', plotFile))
}