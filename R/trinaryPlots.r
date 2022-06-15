#################################################################
#################################################################
#################################################################
#' @title Make Trinary maps from previous workflow runs
#'
#' @description Read in workflow info, make trinary maps and plot vs thresholds
#'
#' @details
#' See Examples.
#'
# #' @param mod.out
# #' @param stats
# #' @param proj.env
# #' @param name
# @keywords
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @export

trinaryMapPlot=function(trinary.rasters,
												plotFile=NULL,
												pres=NULL,
												species,						
												shapesToPlot=NULL,
												openFig=T){

	#  for testing
	#  pres=NULL; shapesToPlot=NULL
	if(!is.null(plotFile))  grDevices::pdf(plotFile,h=7,w=7)

		
		image(trinary.rasters>=0, xlab='',ylab='',col=c('grey80'),xaxt='n', yaxt='n',bty='n',main=species)
		image(trinary.rasters>=1, xlab='',ylab='',col=c('steelblue'),xaxt='n', yaxt='n',bty='n',add=TRUE)
		image(trinary.rasters>=2, xlab='',ylab='',col=c('red1'),xaxt='n', yaxt='n',bty='n',add=TRUE)
		
		if(!is.null(pres)) graphics::points(pres,col='grey90',pch='+',cex=.4)

		if(!is.null(shapesToPlot)) lapply(shapesToPlot,plot,add=T,lwd=.1,border='grey40')

	if(!is.null(plotFile))  grDevices::dev.off(); if(openFig) system(paste0('open ',plotFile))
	
}


#################################################################
#################################################################
#################################################################
#' @title Make Trinary maps from previous workflow runs
#'
#' @description Read in workflow info, make trinary maps and plot vs thresholds
#'
#' @details
#' See Examples.
#'
# #' @param mod.out
# #' @param stats
# #' @param proj.env
# #' @param name
# @keywords
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

#' @export

trinaryROCPlot=function(ROCPlotFile=NULL,
												plotThings,
												out1,
												openFig=T,
												addMapPlot=FALSE){
	#  for testing
	#  plotThings=threshs$plotThings; out1=threshs[[1]]; openFig=T; 												addMapPlot=FALSE

	if(!is.null(ROCPlotFile))  grDevices::pdf(ROCPlotFile,h=1.6*7,w=7.6)
		par(mfrow=c(3,2),mar=c(5,4,2,2))
		plot(plotThings$xx,plotThings$y,type='l',lwd=4,ylim=c(0,1.08), main="ROC",xlab='1-specificity', ylab="sensitivity")
		#abline(h=0,lty=2)
		abline(v=out1$youden.thresh.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		abline(v=out1$lo.thresh.x,col='red1',lty=3,lwd=2)

		abline(0,1,lty=2)
		#points(out1$lo.thresh.x,y.lo,pch=21,bg='red1')
		points(out1$youden.thresh.x,out1$youden.thresh.y,pch=21, bg='darkgoldenrod2')
		points(out1$hi.thresh.x,out1$hi.thresh.y,pch=21,bg='steelblue')
		points(out1$lo.thresh.x,out1$lo.thresh.y,pch=21,bg='red1')
		text(.8,.2,paste0('pAUC =\n',round(plotThings$a.pauc$auc,3)),cex=1.5)
		text(out1$youden.thresh.x+.02,out1$youden.thresh.y-.02, paste0('(',round(out1$youden.thresh.x,2),',', round(out1$youden.thresh.y,2),')'), adj=c(0,0), cex=1,col='darkgoldenrod2') 
		text(out1$hi.thresh.x-.02,out1$hi.thresh.y+.03, paste0('(',round(out1$hi.thresh.x,2),',',round(out1$hi.thresh.y,2),')'),adj=c(1,0), cex=1,col='steelblue')
		text(out1$lo.thresh.x+.02,out1$lo.thresh.y+.03, paste0('(',round(out1$lo.thresh.x,2),',',round(out1$lo.thresh.y,2),')'),adj=c(0,1), cex=1,col='red1')
		
		# ROC' -----------------------------------------
		if(!any(is.finite(plotThings$y.))){ #this catches NAs to
			keep=which(!is.na(plotThings$y.))
			x.tmp=plotThings$xout[keep]; y.tmp=plotThings$y.[keep]
		} 
		p1=plot(x.tmp,y.tmp,type='l',col='black',log='y', main="ROC'", xlab='1-specificity',ylab="sensitivity''",xlim=c(0,1))
		#else {plot(.5,.5,col='white', main="ROC'", xlab='1-specificity',ylab="sensitivity'"); graphics::text(.5,.5,'Derivative contains Inf or NA')}
		abline(h=0,lty=2)
		abline(h=1,lty=2)
		abline(v=out1$youden.thresh.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		abline(v=out1$lo.thresh.x,col='red1',lty=3,lwd=2)

		# ROC'' -----------------------------------------
		if(!any(is.finite(plotThings$y..))){
			keep=which(!is.na(plotThings$y..))
			x.tmp=plotThings$xout[keep]; y.tmp=plotThings$y..[keep]
		} 
		plot(x.tmp,logmod(y.tmp),type='l',col='black', main="logmod(ROC'')", xlim=c(0,1),xlab='1-specificity',ylab="sensitivity")
		#else {
		#	plot(.5,.5,col='white', main="logmod(ROC'')", xlab='1-specificity',ylab="sensitivity'"); graphics::text(.5,.5,'Derivative contains Inf or NA')
		#}
		abline(h=0,lty=2) 
		abline(v=out1$youden.thresh.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		abline(v=out1$lo.thresh.x,col='red1',lty=3,lwd=2)
		
		# inverse ROC -----------------------------------
		if(!any(is.finite(plotThings$y1))){
	 		keep=which(!is.na(plotThings$y1))
			x.tmp=plotThings$xx1[keep]; y.tmp=plotThings$y1[keep]
	 	} 
	 	plot(x.tmp,y.tmp,type='l',lwd=4,ylim=c(0,1.08), main='Inverse ROC',xlab='1-sensitivity',ylab='specificity')
		# 	 	else {
		# 	 		plot(.5,.5,col='white', main="Inverse ROC", xlab='1-sensitivity',ylab='specificity'); graphics::text(.5,.5,'Derivative contains Inf or NA')
		# 	 	}
		abline(h=out1$y.lo.inv,col='red1',lty=3,lwd=2)
		abline(v=out1$x.lo.inv,col='red1',lty=3,lwd=2)
		abline(0,1,lty=2)
		
		# inverse ROC'' -----------------------------------
		if(!any(is.finite(plotThings$y1..))){
			keep=which(!is.na(plotThings$y1..))
			x.tmp=plotThings$x1out[keep]; y.tmp=plotThings$y1..[keep]
		} 
		plot(x.tmp,logmod(y.tmp),type='l', main="logmod[(Inverse ROC)'']",xlab='1-sensitivity',ylab="logmod(specificity' ')",xlim=c(0,1),ylim=range(logmod(plotThings$y1..),na.rm=T))
		# 		else {
		# 		plot(.5,.5,col='white', main="logmod[(Inverse ROC)'']",xlab='1-sensitivity',ylab="logmod(specificity' ')"); graphics::text(.5,.5,'Derivative contains Inf or NA')
		# 		}
		abline(h=y.lo)
		abline(h=0,lty=2)
		abline(v=out1$x.lo.inv,col='red1',lty=3,lwd=2)

		# for third derivative - works but not needed
		# # 		plot(plotThings$xout,logmod(y...),type='l',col='black', log='y',xlim=c(0,1),main="ROC'''")
		# # 		abline(h=0,lty=2)
		# # 		abline(v=x.lo,col='darkgoldenrod2',lty=3,lwd=2)
		# # 		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		# # 		abline(v=x.lo,col='red1',lty=3,lwd=2)

		# curvature; code works, idea doesn't
		# # 		plot(plotThings$xout,c2,type='l',col='black',xlim=c(0,1), main=expression(paste(kappa,"(ROC)")))
		# # 		abline(v=x.lo,col='red1',lty=3,lwd=2)
		# # 		lines(plotThings$xout,c2.1,col='grey30')
		# # 		lines(plotThings$xout,c2.2,col='grey60')
		#abline(v=.y..f(rate,skew),lty=2)
		##abline(v=x.root,col='blue',lty=2)
		#lines(ss$x,ss$y,col='blue')
		
		if(addMapPlot){
			print('not done')
			# have to push all these args, which is annoying
			# trinaryMapPlot(trinary.rasters,plotFile,pres=pres, expertRasterPath=expertRasterPath,expertShpPath=expertShpPath,shapesToPlot=shapesToPlot,openFig=T)
		}

	
	if(!is.null(ROCPlotFile))  grDevices::dev.off()
	if(openFig)	system(paste0('open ', ROCPlotFile))
}