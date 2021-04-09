#################################################################
#################################################################
#################################################################
#' @title Apply trinary map tools to a set of continuous range maps
#'
#' @description Run all aspects of modeling and plotting for a suite of previously run models
#'
#' @details
#' See Examples.
#'
#' @param dirs,
#' @param stats
#' @param env
#' @param pres.f
#' @param plot.f
#' @param map.dir
#' @param plot.expert=TRUE
#' @param openFig=FALSE
#' @param expertRaster.f=NULL

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
trinaryMapWorkflow=function(pres,
													  background,
													  model,
													  modelNames,
													  species,
													  rasterOutputDir=NULL,
													  mapPlotDir=NULL,
													  ROCPlotDir=NULL,
													  expertShpPath=NULL,
													  expertRasterPath=NULL,
													  shapesToPlot=NULL,
													  openFig=T,
													  doMapPlot=TRUE,
													  doROCPlot=TRUE,
													  NATo0=TRUE){
	
	#  for testing
	#  background=bg; doMapPlot=TRUE; doROCPlot=TRUE

	p=raster::extract(model,pres)
	a=raster::extract(model,background)
	if(NATo0){
		a[is.na(a)]=0
		p[is.na(p)]=0
	}
	p <- na.omit(p)
	a <- na.omit(a)
	message(paste0(species,': ',length(p),' presences and ',length(a),' background points used for building trinary maps'))
	ins=rbind(data.frame(Y=1,X=p),data.frame(Y=0,X=a))

	# fit auc curves
	threshs=tryCatch(trinaryROCRoots(ins=ins),error=function(e) e)
	if(class(threshs)=='try-error'){
		message(paste0(species, ": couldn't find roots of the ROC curve; this often happens if you have  very few presence or background points"))
	}
	
	# make maps
	trinary.rasters=trinaryMap(model,
														 threshLo=threshs[[2]]$threshLo,
														 threshHi=threshs[[2]]$threshHi,
														 threshYouden=threshs[[2]]$threshYouden,
														 modelNames=modelNames,
														 species=species,
														 rasterOutputDir=rasterOutputDir,
														 overwrite=TRUE,format="GTiff",
												 		 datatype="INT1U", options=c("COMPRESS=DEFLATE"))
								
	# calculate stats
	range.size=trinaryRangeSize(trinary.rasters)
	
	# plot
	plotFile=paste0(mapPlotDir,modelNames,'.pdf')
	
	if(doMapPlot) trinaryMapPlot(trinary.rasters,plotFile,pres=pres, species=species, expertRasterPath=expertRasterPath,expertShpPath=expertShpPath,shapesToPlot=shapesToPlot,openFig=openFig)
	
	ROCPlotFile=paste0(ROCPlotDir,species,'_',modelNames,'ROC.pdf')
	if(doROCPlot) trinaryROCPlot(ROCPlotFile,plotThings=threshs$plotThings,out1=threshs[[1]],openFig=openFig)
	
	return(list(threshs=threshs,trinary.rasters=trinary.rasters))
	
}