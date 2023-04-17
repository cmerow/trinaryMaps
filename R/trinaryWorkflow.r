#################################################################
#################################################################
#################################################################
#' Generate trinary maps from species presence, background and model
#'
#' This function takes in the species presence and background data, along with a model, and generates trinary maps. It also calculates the range size and plots the maps and ROC curves.
#'
#' @param pres A raster object containing the species presence data
#' @param background A raster object containing the background data
#' @param model A raster object containing the model data
#' @param main An optional title for the main plot
#' @param rasterOutputPath A character string specifying the path to save the trinary rasters
#' @param mapPlotPath A character string specifying the path to save the trinary map plot
#' @param ROCPlotPath A character string specifying the path to save the ROC curve plot
#' @param maxTPQuantile A numeric value between 0 and 1 specifying the maximum true positive quantile to use for threshold selection
#' @param shapesToPlot A list containing one or more SpatialPolygonsDataFrame objects to overlay on the map plot
#' @param openFig A logical indicating whether the generated plots should be displayed in the R console
#' @param NATo0 A logical indicating whether to turn NA values into zeros in the presence and background data
#'
#' @return A list containing the trinary thresholds and trinary rasters
#'
#' @export
trinaryMapWorkflow=function(pres,
													  background,
													  model,
													  main=NULL,
													  #modelNames,
													  #species,
													  #rasterOutputDir=NULL,
													  #mapPlotDir=NULL,
													  #ROCPlotDir=NULL,
													  rasterOutputPath=NULL,
													  mapPlotPath=NULL,
													  ROCPlotPath=NULL,
													  maxTPQuantile,
													  shapesToPlot=NULL,
													  openFig=T,
													  NATo0=TRUE){
	
	#  for testing
	#  background=bg1; ROCPlotDir=mapPlotDir

	p=raster::extract(model,pres)
	a=raster::extract(model,background)
	if(NATo0){ #optionally turn NA into zeros so they contribute to getting absences right. i chose to do this because bg from neightboring ecoregions are considered absences for this evaluation
		a[is.na(a)]=0
		p[is.na(p)]=0
	}
	p <- stats::na.omit(p)
	a <- stats::na.omit(a)
	message(paste0(length(p),' presences and ',length(a),' background points used for building trinary maps'))
	ins=rbind(data.frame(Y=1,X=p),data.frame(Y=0,X=a))

	#== fit auc curves
	threshs=tryCatch(trinaryROCRoots(ins=ins),error=function(e) e)
	if(class(threshs)=='try-error'){
		message(paste0("Couldn't find roots of the ROC curve; this often happens if you have  very few presence or background points. So you're not getting any trinary maps."))
		return(list(threshs=NULL,trinary.rasters=NULL))
	}
	
	#== make maps
	trinary.rasters=trinaryMap(model,
														 threshLo=threshs[[1]]$threshLo,
														 threshHi=threshs[[1]]$threshHi,
														 rasterOutputPath=rasterOutputPath,
														 overwrite=TRUE,format="GTiff",
												 		 datatype="INT1U", options=c("COMPRESS=DEFLATE"))
								
	#== calculate stats
	range.size=trinaryRangeSize(trinary.rasters,model>threshs[[1]]$threshYouden)
	
	#== plot
	if(!is.null(mapPlotPath)) trinaryMapPlot(trinaryRaster=trinary.rasters,plotFile=mapPlotPath,  pres=pres, main=main, shapesToPlot=shapesToPlot,openFig=openFig)
	
	if(!is.null(ROCPlotPath)) trinaryROCPlot(trinaryPlotThings=threshs$plotThings,trinaryDF=threshs$trinaryDF,plotFile=ROCPlotPath,openFig=openFig)
	
	return(list(threshs=threshs,trinary.rasters=trinary.rasters))
	
}