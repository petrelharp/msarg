#' Check a demography object to see if final states are communicating.
#'
#' @param dem The demography object.
#' @param ga The time step to check.
#' @param nstarts The number of randomly chosen starting locations to test.
#' @param niter The number of forward iterations of the migration matrix.
#' @export
check_demography <- function (dem,
                              ga=dem[[length(dem)]],
                              nstarts=10,
                              niter=1000) {
    # check if the longest-ago setup communicates between states
    # by spreading out mass from from random initial states
    M <- ga@M
    M <- M/(2*max(Matrix::rowSums(M))) + Matrix::Diagonal(nrow(M),1-Matrix::rowSums(M))
    xx <- do.call( rbind, lapply( 1:nstarts, function (k) {
                        ifelse( 1:nrow(M)==sample.int(nrow(M),1), 1, 0 ) } ) )
    for (k in 1:niter) {
        xx <- xx %*% M
    }
    return( ! any( Matrix::colSums(xx) == 0 ) )
}


