#################################################################
#################################################################
#################################################################
#' @export
trinaryMapPlot=function(plotF=NULL,trinary.rasters,pres,expertRaster.f=NULL,expertShp=NULL,shapesToPlot=NULL,openFig=T){

	if(!is.null(plotF)){
		pdf(plotf,h=7,w=7)
	}

		image(trinary.rasters, xlab='',ylab='',col=c('grey80','steelblue','red1'),xaxt='n', yaxt='n',bty='n',main=species)
		if(!is.null(expertRaster.f)){
			true=raster(expertRaster.f)
			image(true,add=T,col=c(rgb(0,0,0,0),rgb(0,0,0,.5)))
		}
		points(pres,col='grey90',pch='+',cex=.4)
		if(!is.na(expertShp)) {
			 tmp=readShapePoly(expertShp)
			 projection(tmp)=projection(trinary.rasters)
			 plot(tmp,add=T,border="black",col=grey(.6,.3),lwd=.9)
		}
		if(!is.null(shapesToPlot)){
			lapply(shapesToPlot,plot,add=T,lwd=.1,border='grey40')
		}

	if(!is.null(plotF)){
		dev.off()
		if(openFig) system(paste0('open ',plotf))
	}
}


#======================================================================
#======================================================================

#' @export

trinaryROCPlot=function(plot.f=NULL,plotThings,out1){
	if(!is.null(plot.f)) pdf(plot.f,h=1.6*7,w=7.6)
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
		text(.8,.2,paste0('pAUC =\n',round(a.pauc$auc,3)),cex=1.5)
		text(out1$youden.thresh.x+.02,out1$youden.thresh.y-.02, paste0('(',round(out1$youden.thresh.x,2),',', round(out1$youden.thresh.y,2),')'), adj=c(0,0), cex=1,col='darkgoldenrod2') 
		text(out1$hi.thresh.x-.02,out1$hi.thresh.y+.03, paste0('(',round(out1$hi.thresh.x,2),',',round(out1$hi.thresh.y,2),')'),adj=c(1,0), cex=1,col='steelblue')
		text(out1$lo.thresh.x+.02,out1$lo.thresh.y+.03, paste0('(',round(out1$lo.thresh.x,2),',',round(out1$lo.thresh.y,2),')'),adj=c(0,1), cex=1,col='red1')
		
		plot(plotThings$xout,plotThings$y.,type='l',col='black',log='y',main="ROC'", xlab='1-specificity',ylab="sensitivity'")
		abline(h=0,lty=2)
		abline(h=1,lty=2)
		abline(v=out1$youden.thresh.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		abline(v=out1$lo.thresh.x,col='red1',lty=3,lwd=2)

		plot(plotThings$xout,logmod(plotThings$y..),type='l',col='black', main="logmod(ROC'')", xlim=c(0,1),xlab='1-specificity',ylab="sensitivity' '")
		abline(h=0,lty=2) 
		abline(v=out1$youden.thresh.x,col='darkgoldenrod2',lty=3,lwd=2)
		abline(v=out1$hi.thresh.x,col='steelblue',lty=3,lwd=2)
		abline(v=out1$lo.thresh.x,col='red1',lty=3,lwd=2)
		
	 	plot(plotThings$xx1,plotThings$y1,type='l',lwd=4,ylim=c(0,1.08),main='Inverse ROC',xlab='1-sensitivity',ylab='specificity')
		abline(h=out1$y.lo.inv,col='red1',lty=3,lwd=2)
		abline(v=out1$x.lo.inv,col='red1',lty=3,lwd=2)
		abline(0,1,lty=2)
		
		plot(plotThings$x1out,logmod(plotThings$y1..),type='l',main="logmod[(Inverse ROC)'']",xlab='1-sensitivity',ylab="logmod(specificity' ')",xlim=c(0,1),ylim=range(logmod(plotThings$y1..),na.rm=T))
		#abline(h=y.lo)
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
		if(!is.null(trinaryPlotArgs)){
			range.size=trinaryMap(dirs=dirs, modelPath=trinaryPlotArgs$modelPath, threshLo=plotThings$threshLo,threshHi=plotThings$threshHi,threshYouden=plotThings$threshYouden,expertShp=trinaryPlotArgs$expertShp, shapesToPlot=trinaryPlotArgs$shapesToPlot, pres=trinaryPlotArgs$pres,map.dir=plotThings$map.dir,expertRaster.f=trinaryPlotArgs$expertRaster.f)$range.size
		} else {
			range.size=NULL
		}
	
	if(!is.null(plot.f)) dev.off()
	if(openFig)	system(paste0('open ', plot.f))
	range.size
}