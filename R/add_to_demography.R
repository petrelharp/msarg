
#' @export
add_to_demography <- function (dem,dt,fn=identity,pop,...,tnew=dt+if(length(dem@t)>0){dem@t[length(dem@t)]}else{0}) {
    # add to a demography by applying a modifier function fn 
    #  to another population model,
    #  by default the previous state
    if (inherits(dem,"popArray")) {
        dem <- new("demography",popStates=list(dem),t=numeric(0))
    }
    nt <- length(dem@popStates)
    if (missing(pop)) { pop <- dem@popStates[[nt]] }
    if (is.numeric(pop)) { pop <- dem@popStates[[pop]] }
    dem@popStates <- c( dem@popStates, fn(pop,...) )
    stopifnot( tnew >= max(c(0,dem@t)) )
    dem@t <- c(dem@t,tnew)
    return(dem)
}

