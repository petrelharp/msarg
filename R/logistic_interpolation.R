#' Modify a Demography by Interpolating a Logistic Expansion
#'
#' Does logistic growth to interpolate between population states.
#'
#' @param dem The demography object.
#' @param t.end The time ago to end the interpolation.
#' @param t.begin The time ago to begin the interpolation (should be *larger* than t.end).
#' @param nsteps The number of time steps to add to the demography.
#' @param speed The speed of the expanding wave, in demes per unit time.
#' @param width The width of the expanding wave, in demes.
#' @param sigma The smoothing factor in the logistic wave: a dispersal distance per unit time.
#' @param r The logistic rate of change.
#' @param dt.per.step Number of 'invisible' steps between each recorded step (should be fairly large).
#' @param max.step Largest (invisible) step size, as a proportion, for population sizes.
#' @export
logistic_interpolation <- function (dem, 
                                    t.end, 
                                    t.begin, 
                                    nsteps, 
                                    speed, 
                                    width, 
                                    sigma=sqrt(speed*width/sqrt(2)), 
                                    r=speed/(width*sqrt(2)), 
                                    dt.per.step=100,
                                    max.step=0.2) {
    # Do logistic growth to interpolate between population states.
    #  does dt.per.step 'invisible' steps -- it seems essential this is fairly large.
    stopifnot( ( t.end < t.begin ) && ( t.begin <= max(dem@t) ) )
    #   This only works if both t.end and t.begin are in the same interval,
    #     with population states specified on either end.
    ga.interval <- findInterval( (t.end+t.begin)/2, c(0,dem@t) )
    stopifnot( ( c(0,dem@t)[ga.interval] <= t.end ) && ( t.begin <= c(0,dem@t)[ga.interval+1] ) ) 
    # we will add new population states at these times:
    add.t <- seq(t.begin,t.end,length.out=nsteps+1)
    new.t <- sort(unique(c(dem@t,add.t)))
    # full list of demographies will go here
    new.dem <- vector(length(new.t),mode="list")
    # populate with existing states
    new.dem[ match(c(0,dem@t),c(0,new.t)) ] <- dem@popStates
    # ending state, determines carrying capacities
    end.N <- dem[[ga.interval]]@N
    # beginning state
    begin.k <- match(t.begin,c(0,new.t))
    new.dem[[ begin.k ]] <- next.ga <- dem[[ga.interval+1]]
    zeros <- ( (next.ga@N==0) & (end.N==0) )
    # time step length
    dt <- (t.begin-t.end)/(nsteps*dt.per.step)
    # Laplacian operator -- should really use *migration matrix* for this?
    #  do some gymnastics to get migration matrix if sigma > 1
    nsig <- max(0,ceiling(log2(dt*sigma^2)))
    sig2 <- (dt*sigma^2) / 2^nsig
    adj <- (sig2/4) * grid.adjacency(nrow(dem),ncol(dem),diag=FALSE)
    diag(adj) <- 1-Matrix::rowSums(adj)
    for (k in seq_len(1+nsig)[-1]) { adj <- adj%*%adj }
    for (k in 1:nsteps) {
        for (j in 1:dt.per.step) {
            ## discrete logistic, but with bounded steps
            lstep <- pmax(-max.step, pmin(max.step, dt * r * ( end.N - next.ga@N ) ) )
            next.ga@N <- as.vector( adj %*% ( next.ga@N * ( 1 + lstep ) ) )
            next.ga@N[zeros] <- 0.0
        }
        new.dem[[ begin.k - k ]] <- next.ga
    }
    return( 
        new( "ms_demog",
            popStates=new.dem,
            t=new.t
        ) )
}


