setClass("NicheFun", representation(shape = "character", location = "numeric", breadth = "numeric", 
								 scale = "numeric", nr = "integer"), validity = check_niche)

check_niche_fun = function(x) {
	errors = character()
	if(length(shape) != 1 || !shape %in% c("gaussian", "constant", "inverse_gaussian"))
		errors = c(errors, "Shape must be gaussian, inverse_gaussian or constant")
	if(shape == "gaussian" | shape == "inverse_gaussian") {
		if(length(location) != nr)
			errors = c(errors, "Must have one location per resource dimension")
		if(any(breadth < 0))
			errors = c(errors, "Niche breadth must be positive")
		if(nr == 1) {
			if(length(br != 1))
				errors = c(errors, "Single-dimensional niche functions must have a single niche breadth")
		} else {
			if(!is(breadth, "matrix") && nrow(breadth) != nr && ncol(breadth) != nr)
				errors = c(errors, "Multi-dimensional niche breadth must be a matrix with nrow = ncol = nr")
		}
		if(length(scale) != 1)
			errors = c(errors, "Scale must be a single value")
	}
	if(length(errors) == 0) TRUE else errors
}
