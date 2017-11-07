#' @export
logmod=function(x){sign(x)*log(abs(x))}


#' @export
# .curv1=function(y.,y..,y...){
# 	abs(y.*y..)/abs(y...)
# }

#' @export
.curv2=function(y.,y..){
	abs(y..)/( (1+(y.^2))^1.5 )
}
