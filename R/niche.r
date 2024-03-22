#' @keywords internal
.check_niche = function(object) {
	errors = character()
	if(length(object@scale) != 1)
		errors = c(errors, "Scale must be a single value")
	if(length(errors) == 0) TRUE else errors
}

#' @keywords internal
.check_gniche = function(object) {
	errors = character()
	if(length(object@location) != object@nr)
		errors = c(errors, "Must have one location per resource dimension")
	if(any(object@breadth < 0))
		errors = c(errors, "Niche breadth must be positive")
	if(nrow(object@breadth) != object@nr && ncol(object@breadth) != object@nr)
		errors = c(errors, "Niche breadth must be a matrix with nrow == ncol == nr")
	if(length(errors) == 0) TRUE else errors
}

# #' @keywords internal
# .check_igniche = function(object) {
# 	errors = character()
# 	if(length(object@displacement) != 1)
# 		errors = c(errors, "Displacement must be a single value")
# 	if(length(errors) == 0) TRUE else errors
# }

setClass("NicheFun", representation(scale = "numeric", nr = "integer"), validity = .check_niche)
setClass("ConstantNicheFun", contains = "NicheFun")
setClass("GaussianNicheFun", representation(location = "numeric", breadth = "matrix"), 
		 contains = "NicheFun", validity = .check_gniche)
## To implement: right now we do this via ce_funs. These should move over to be S4 methods
# setClass("InverseGaussianNicheFun", representation(displacement = "numeric"), 
# 		 contains = "GaussianNicheFun", validity = .check_igniche)

#' Get niche Parameters
#' @name get_pars
#' @keywords internal 
setGeneric("get_pars", function(object) {
	standardGeneric("get_pars")
})

#' @rdname get_pars
#' @keywords internal 
setMethod("get_pars", signature("NicheFun"), function(object) list(scale = object@scale))

#' @rdname get_pars
#' @keywords internal 
setMethod("get_pars", signature("ConstantNicheFun"), function(object) {
	c(list(nr = object@nr), callNextMethod())
})

#' @rdname get_pars
#' @keywords internal 
setMethod("get_pars", signature("GaussianNicheFun"), function(object) {
	c(list(location = object@location, breadth = object@breadth), callNextMethod())
})
# setMethod("get_pars", signature("InverseGaussianNicheFun"), function(object) {
# 	c(list(displacement = object@displacement), callNextMethod())
# })
