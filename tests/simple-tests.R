#!/usr/bin/Rscript

library(msarg)

##
# a 2x3 grid, sampling at three corners:
#
#  1 . .
#  2 . 3
#
# and the demes are ordered like:
#
#  1 3 5
#  2 4 6
#
# and having sizes
#
#  2 4 6
#  1 3 5

deme.sizes <- c(2,1,4,3,6,5)
ga <- grid_array(nlayers=1,nrow=2,ncol=3,N=deme.sizes,mig.rate=10)
sample.config <- sort_sample_config( cbind(
        row=c(1,2,2),
        col=c(1,1,3),
        layer=1,
        n=c(1,2,3)
    ) )
dem <- ms_demog( ga )

# so sample vector should be 1,2,0,0,0,3
nsvec <- nsample_vector(dem[[1]],sample.config)
stopifnot( all.equal( nsvec, c(1,2,0,0,0,3) ) )

# and the command-line for ms is
ms.arg <- msarg(dem,nsvec,scale.migration=FALSE)

# sampling vector
#   -I 6 1 2 0 0 0 3
stopifnot( all( ms.arg[["-I"]] == c(length(nsvec),nsvec) ) )

# relative population sizes
#   -n 1 1 -n 2 1 -n 3 1 -n 4 1 -n 5 1 -n 6 1
ms.sizes <- do.call(rbind,ms.arg[names(ms.arg)=="-n"])
stopifnot( all( ms.sizes[ms.sizes[,1],2] == deme.sizes ) )

# migration rates:
#   -m 2 1 10 -m 3 1 10 -m 1 2 10 -m 4 2 10 -m 1 3 10 -m 4 3 10 -m 5 3 10 -m 2 4 10 -m 3 4 10 -m 6 4 10 -m 3 5 10 -m 6 5 10 -m 4 6 10 -m 5 6 10
#    1 3 5
#    2 4 6
true.m <- matrix( c(
        1, 2,
        1, 3,
        2, 1,
        2, 4,
        3, 1,
        3, 4,
        3, 5,
        4, 2,
        4, 3,
        4, 6,
        5, 3,
        5, 6,
        6, 4,
        6, 5
    ), ncol=2, byrow=TRUE )
# true.m <- cbind( true.m, 10 * deme.sizes[true.m[,1]]/deme.sizes[true.m[,2]] )
true.m <- cbind( true.m, 10 )
true.m <- true.m[ order(true.m[,1],true.m[,2]), ]
ms.m <- do.call(rbind,ms.arg[names(ms.arg)=="-m"])
ms.m <- ms.m[ order(ms.m[,1],ms.m[,2]), ]
stopifnot( all.equal( ms.m, true.m, check.attributes=FALSE ) )

##
# check tree things
msout <- run_ms( dem, sample.config, tofile=FALSE, trees=TRUE, seeds="1 2 3" )
trees <- trees_from_ms(msout)
dists <- tree_dists( trees, sample.config )
stopifnot( all( colnames(dists[[1]]) == as.character(1:6) ) && all( rownames(dists[[1]]) == as.character(1:6) ) )


##
# a 5x4 grid, sampling at three corners:
#
#  1 . . 3
#  . . 2 .
#  . . . .
#  . . . .
#  5 . . .
#
# and the demes are ordered like:
#
#  1 6  11 16 
#  2 7  12 17 
#  3 8  13 18 
#  4 9  14 19 
#  5 10 15 20 
#

ga <- grid_array(nlayers=1,nrow=5,ncol=4,N=1,mig.rate=10)
sample.config <- sort_sample_config( cbind(
        row=c(1,5,2,1),
        col=c(1,1,3,4),
        layer=1,
        n=c(1,5,2,3)
    ) )
dem <- ms_demog( ga )

# so sample vector should be 1,2,0,0,0,3
nsvec <- nsample_vector(dem[[1]],sample.config)
stopifnot( all.equal( nsvec, c(1,0,0,0,5,0,0,0,0,0,0,2,0,0,0,3,0,0,0,0) ) )

