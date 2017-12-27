
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
logmod=function(x){sign(x)*log(abs(x))}


# @export
# .curv1=function(y.,y..,y...){
# 	abs(y.*y..)/abs(y...)
# }

#' @export
.curv2=function(y.,y..){
	abs(y..)/( (1+(y.^2))^1.5 )
}
