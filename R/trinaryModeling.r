#################################################################
#################################################################
#################################################################
#' @title Find upper/lower thresholds from ROC curves
#'
#' @description Fit a smoothed ROC curve, find bounds for threshold and report partial AUC statistics
#'
#' @details
#' See Examples.
#'
#' @param ins
# @param plot.f
#' @param max.sens=0.95
# @param doPlot=FALSE
# @param trinaryPlotArgs
# @param openFig=FALSE
# @param mapOutputDir=NULL
# @param ... options to be passed to `pROC::smooth`
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
trinaryROCRoots=function(ins,
 												 max.sens=0.95,
 												 smoothMethod='binormal',...){
 												 
out=try({
	
	aa=pROC::roc(ins[,1], ins[,2],quiet=T)
	a=try(pROC::smooth(aa,method=smoothMethod, n=1024,...),silent=TRUE)
	# i used this default because its the pROC package defualt so i assumed it was the best. if it breaks, try the next one
	if(class(a)=='try-error') { 
		a=pROC::smooth(aa,method='density',...)
		message("Used method=density for ROC smoothing because your selected method (set by argument smoothMethod) broke. If you're unhappy about this, see other options for methods in ?pROC::smooth." )
	}
	xx=rev(1-a$specificities)
	xx.=.middle_pts(xx)
	xx..=.middle_pts(.middle_pts(xx))
	xx...=.middle_pts(.middle_pts(.middle_pts(xx)))

	xout=seq(0,1,length=200)
	y=rev(a$sensitivities)
	y.r= .deriv(xx, y)
	y..r <- .deriv(xx., y.r)
	y...r <- .deriv(xx.., y..r)
	y....r <- .deriv(xx..., y...r)
	
	#== make the derivatives real functions so they can be evaluated at the x points (e.g. for curvature)
	y.=suppressWarnings(approx(xx.,y.r,xout=xout,method='linear')$y)
	y..=approx(xx..,y..r,xout=xout,method='linear')$y
	y...=approx(xx...,y...r,xout=xout)$y
	y....=approx(.middle_pts(xx...),y....r,xout=xout)$y
	
	#== youden index 
	youden=a$specificities+a$sensitivities-1
	best.youden=which.max(youden)
	y.youden=a$sensitivities[best.youden]
	x.youden=(1-a$specificities)[best.youden]
	
	#== remove NAs 
	keep=complete.cases(y..r)
	y..1=y..r[keep]
	xx..1=xx..[keep]

	#== get locations of sign changes of the logmod of the second derivative
	(switches=cumsum(rle(sign(logmod(y..1)))[[1]]))
	#-- remove  nans
	switches[which(switches==length(y..1))]=NA
	switches=na.omit(switches)
	
		#-- since y is evaluated at the midpoint of the xs get the midpoint...
	if(length(switches)>0 & any(switches>best.youden)){
		x.as=(xx..1[switches]+xx..1[switches+1])/2
		#### x.as=xout[switches]
		x.as=x.as[which(x.as>x.youden)[1]] # this finds the midpoint, need to find the left interval
		
		x.ind=findInterval(x.as,xx)
		y.as=y[x.ind]

	} else {
		#-- if no asymptote reached, use max sens
		y.as=max.sens
		not.root.index=findInterval(max.sens,y)
		x.as=xx[not.root.index]
	}

	#== prep for COR (inverse ROC) to find asymptote
	y1=1-xx
	xx1=1-y
	xx1.=.middle_pts(xx1)
	xx1..=.middle_pts(.middle_pts(xx1))
	xx1...=.middle_pts(.middle_pts(.middle_pts(xx1)))
	x1out=seq(0,1,length=200)
	y1.r= .deriv(xx1, y1)
	y1..r <- .deriv(xx1., y1.r)
	y1...r <- .deriv(xx1.., y1..r)
	y1....r <- .deriv(xx1..., y1...r)
	y1.=suppressWarnings(approx(xx1.,y1.r,xout=xout)$y)
	y1..=approx(xx1..,y1..r,xout=xout)$y
	y1...=approx(xx1...,y1...r,xout=xout)$y
	y1....=approx(.middle_pts(xx1...),y1....r,xout=xout)$y
	####keep=complete.cases(y1..)
	####y1..1=y1..[keep]
	keep=complete.cases(y1..r)
	y1..1=y1..r[keep]
	(switches=cumsum(rle(sign(logmod(y1..1)))[[1]]))
	#-- remove  nans
	switches[which(switches==length(y1..1))]=NA
	switches=na.omit(switches)
	
		#-- since y is evaluated at the midpoint of the xs get the midpoint...
	#if(length(switches)>0 & any(switches>best.youden)){
	x.lo.inv=min((xx1..[switches]+xx1..[switches+1])/2)
	#### x.lo=max(xout[switches]) #- not sure this will generally work 
	#x.lo=x.lo[which(x.lo<y.lo)[1]] # just hoping this is ok
	
	x.ind=findInterval(x.lo.inv,rev(xx1))
	y.lo.inv=rev(y1)[x.ind]

	#} else {
	# 		#-- if no asymptote reached, use max sens
	# commented because hopefully not an issue
	# 		y.as=max.sens
	# 		not.root.index=findInterval(max.sens,y)
	# 		x.as=xx[not.root.index]
		#}
		# plot(xx1,y1,type='l')
	# 	abline(h=y.lo)
	# 	abline(v=x.lo)
	# 	plot(x1out,y1.,type='l',log='y')
	# 	plot(x1out,y1..,type='l',log='y')
	
	x.lo=1-y.lo.inv
	y.lo=1-x.lo.inv

	a.pauc=try(pROC::roc(ins[,1], ins[,2],auc=T,partial.auc=1-c(x.lo,x.as), partial.auc.focus='specificity', partial.auc.correct=T,smooth=TRUE,quiet=T),silent=T)
	if(is.na(a.pauc)) a.pauc=try(pROC::roc(ins[,1], ins[,2],auc=T,partial.auc=1-c(x.lo,x.as), partial.auc.focus='specificity', partial.auc.correct=T,smooth.method='density',quiet=T),silent=T)
	if(class(a.pauc)=='try-error') a.pauc=list(auc=NA)
	
	#== find thresholds 
	#-- smoothing the auc, which was needed for derivativies, doesn't give you thresholds associated with the prediction (it just smooths the ROC curve, and there are no underlying threshold values associated with this smooth) so i ask for the data threshold associated with the smoothed curve. It can happen that when there are few points, the same threshold is associated with different hi, youden, low estimates that came from smoothing. This  means that all the maps are the same.
	#-- cm: 9/3/21 i canned this approach and used the smoothed curves because in cases where there are big jumps between data points or AUC ~1, you don't get any differences in the estimated threoushold 
	a.rough=pROC::roc(ins[,1], ins[,2],quiet=T)
	threshLo=rev(a.rough$thresholds)[findInterval(x.lo, rev(1-a.rough$specificities))]
	threshYouden=rev(a.rough$thresholds)[findInterval(x.youden, rev(1-a.rough$specificities))]
	threshHi=rev(a.rough$thresholds)[findInterval(x.as, rev(1-a.rough$specificities))]
	
	#== prep outputs
	out1=data.frame(lo.thresh.x=x.lo,           lo.thresh.y=y.lo, 			
									youden.thresh.x=x.youden,   youden.thresh.y=y.youden,
									hi.thresh.x=x.as,           hi.thresh.y=y.as,
									y.lo.inv=y.lo.inv,          x.lo.inv=x.lo.inv,
									trinary.pauc=as.numeric(a.pauc$auc))
	plotThings=list(xx=xx,y=y,y.=y.,y..=y..,xx1=xx1,y1=y1,x1out=x1out,y1..=y1.., a.pauc=a.pauc,xout=xout,x1out=x1out,threshLo=threshLo,threshHi=threshHi,threshYouden=threshYouden)
	
	list(out1=out1,plotThings=plotThings)
	
	}) 
	return(out)
}



#################################################################
#################################################################
#################################################################
#' @title Make trinary maps 
#'
#' @description Use previously calculated thresholds to make trinary maps
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

# stats to grab model paths (could move that outside as an input arf
# inputPaths=stats[!is.na(stats$meanModelPath),] (could use modelPath already)
# dirs to get species name and define output directory
# outputPaths=paste0(
# !!!!!! shoulds split this into making the map and plotting it.

# should add the option to use a map in memory
trinaryMap=function(model,
										threshLo,
										threshHi,
										threshYouden,
										modelNames,
										species,
										rasterOutputDir=NULL,
										...
										){
										
	out=try({	
	 
	 #== set up inputs for workflow use or independent use
	# 	 if(is.null(modelPath)){
	# 	 	 species=basename(dirs$sp.pred.path)
	# 	 	 stats=stats[!is.na(stats$meanModelPath),]
	# 	 	 preds.r=cmsdm::unrtrans(raster::stack(stats$meanModelPath))	 
	# 	 	 modelName=stats$modelName[x]
	# 	 } else {
	# 	 	 species=basename(dirname(modelPath))
	# 	 	 preds.r=cmsdm::unrtrans(raster::raster(modelPath))
	# 	 	 modelName=tools::file_path_sans_ext(basename(modelPath))
	# 	 }
	preds.r=model 
	#== make trinary maps
	trinary.rasters=stack(lapply(1:nlayers(preds.r),
												function(x) { 
													#stack(lapply(1:ncol(thresh[x,]),
													 #function(y) {
													 out1=preds.r[[x]]>threshLo
													 out2=preds.r[[x]]>threshHi
													 out3=out1+out2
													 names(out3)=paste0(modelNames[x],'_', round(threshLo,2),'_',round(threshHi,2))
													 out3
											}))#}))	 
									 
	 #plot(hist(log(values(preds.r))))
	 #abline(v=log(threshLo)); abline(v=log(threshHi))
	  
		# just make this to get the range size
  youden.binary.tmp.raster=preds.r>threshYouden
											
	 #== write results
	if(!is.null(rasterOutputDir)){
		 lapply(1:nlayers(trinary.rasters),function(x){
			 if(is.null(rasterOutputDir)){
			 	trinaryDirSp=paste0(dirs$trinaryDir,'/',species)
			 	if(!file.exists(trinaryDirSp)) dir.create(trinaryDirSp)
			 	trinaryModel=paste0(trinaryDirSp,'/', names(trinary.rasters)[x] ,".tif")
			 } else {
			 		trinaryModel=paste0(rasterOutputDir,'/',species,'_', names(trinary.rasters)[x] ,".tif")
			 }
			 
			 rf <- suppressWarnings(writeRaster(trinary.rasters[[x]], 
												 filename=trinaryModel, ...))
		 })
	  } # end if is.null(rasterF) 	

	  return(trinary.rasters)
	})
	return(out)
}

#################################################################
#################################################################
#################################################################
#' @title Calculate upper and lower limits of range size
#'
#' @description Size limits based on trinary thresholds
#'
#' @details
#' See Examples.
#'
#' @param trinary.rasters,
#' @param youden.binary.tmp.raster
# @keywords
#'
# @examples
#'
#'
#' @return a data.frame with the lower, upper, and optimal (Youden threshold) range size
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

#' @export
trinaryRangeSize=function(trinary.rasters,
                          youden.binary.tmp.raster=NULL){
  
  cell.size=prod(res(trinary.rasters)/1e3)
	if(!is.null(youden.binary.tmp.raster)){ range.size.youden.km2=sum(values(youden.binary.tmp.raster)>0,na.rm=T) 
	} else {
		range.size.youden.km2=NA
	}
	
	range.size= cell.size * data.frame( 
				range.size.lo.km2=sum(values( trinary.rasters)>1,na.rm=T),
	 			range.size.hi.km2=sum(values( trinary.rasters)>0,na.rm=T),
	 			range.size.youden.km2=range.size.youden.km2)

	range.size
}
	


