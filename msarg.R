
require(Matrix)

##
# A demographic model consists of:
#  - a set of populations (demes)
#  - their sizes
#  - and their pairwise migration rates
# all possibly varying in time, and so hence also:
#  - an identification of populations between time points

# An entire demographic model strings together states:
setClass("demography", representation(
        popStates="list",       # a list of popArrays
        t="numeric"             # vector of lengths of times ago of transitions between the states
    ) )

setMethod("length", signature=c(x="demography"), definition=function (x) { length(x@t) } )
setMethod("[[", signature=c(x="demography",i="ANY"), definition=function (x,i) { x@popStates[[i]] } )
setMethod("[[<-", signature=c(x="demography",i="ANY",value="ANY"), definition=function (x,i,value) { x@popStates[[i]]<-value; return(x) } )

add_to_demography <- function (dem,dt,fn,...,tnew=if(length(dem@t)>0){dem@t[length(dem@t)]}else{0}+dt) {
    # add to a demography by applying a modifier function fn to the previous state
    if (inherits(dem,"popArray")) {
        dem <- new("demography",popStates=list(dem),t=numeric(0))
    }
    nt <- length(dem@popStates)
    dem@popStates <- c( dem@popStates, fn(dem@popStates[[nt]],...) )
    dem@t <- c(dem@t,tnew)
    return(dem)
}

pop_to_dem <- function (ga,t) { new("demography",popStates=list(ga),t=numeric(0)) }

setClass("msarg", contains="namedList")
setMethod("print", signature=c(x="msarg"), definition=function (x,collapse="\n") {
        cat( paste( paste(names(x),sapply(lapply(x,unlist),paste,collapse=' ')), collapse=collapse ), "\n" )
    } )

setGeneric("msarg", function(x,...) { standardGeneric("msarg") })
setMethod("msarg", signature=c(x="demography"), definition=function (x,theta,...) {
        ## FROM ms:
        # usage: ms nsam howmany 
        #   Options: 
        # 	 -t theta   (this option and/or the next must be used. Theta = 4*N0*u )
        # 	 -s segsites   ( fixed number of segregating sites)
        # 	 -T          (Output gene tree.)
        # 	 -F minfreq     Output only sites with freq of minor allele >= minfreq.
        # 	 -r rho nsites     (rho here is 4Nc)
        # 		 -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) 
        # 			 if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.
        # 	 -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t
        # 	 -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) 
        # 		 -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)
        # 		 -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)
        # 		 -n i size_i   (popi has size set to size_i*N0 
        # 		 -g i alpha_i  (If used must appear after -M option.)
        # 	   The following options modify parameters at the time 't' specified as the first argument:
        # 	 -eG t alpha  (Modify growth rate of all pop's.)
        # 	 -eg t i alpha_i  (Modify growth rate of pop i.) 
        # 	 -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)
        # 	 -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )
        # 	 -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)
        # 	 -eN t size  (Modify pop sizes. New sizes = size*N0 ) 
        # 	 -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)
        # 	 -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.
        # 		 proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.
        # 		 Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.
        # 	 -ej t i j   ( Join lineages in pop i and pop j into pop j
        # 		  size, alpha and M are unchanged.
        # 	  -f filename     ( Read command line arguments from file filename.)
        # 	  -p n ( Specifies the precision of the position output.  n is the number of digits after the decimal.)
        arglist <- c( list(
                        "-t"=theta,
                        "-I"=c( prod(x[[1]]@npop), x[[1]]@N )
                        ),
                    msarg_G(x[[1]]@G),
                    msarg_M(x[[1]]@M)
                )
        for (k in seq_along(x@t)) {
            t <- x@t[k]
            arglist <- c( 
                    arglist,
                    msarg_N(x[[k+1]]@N,t,old.N=x[[k]]@N),
                    msarg_M(x[[k+1]]@M,t),
                    msarg_G(x[[k+1]]@G,t,old.G=x[[k]]@G)
                )
        }
        return(new("msarg",.Data=arglist,names=names(arglist)))
    } )

msarg_N <- function (N,t=numeric(0),old.N) {
    # Ns that haven't changed don't matter
    # so only output the relevant -eg commands
    dothese <- if (missing(old.N)) { seq_along(N) } else { which(N!=old.N) }
    Narg <- lapply( seq_along(dothese), function (k) {
            c( t, dothese[k], N[k] )
        } )
    names( Narg ) <- rep(if (length(t)>0){"-en"}else{"-n"},length(Narg))
    return(Narg)
}

msarg_G <- function (G,t=numeric(0),old.G) {
    # Gs that haven't changed don't matter
    # so only output the relevant -eg commands
    dothese <- if (missing(old.G)) { seq_along(G) } else { which(G!=old.G) }
    Garg <- lapply( seq_along(dothese), function (k) {
            c( t, dothese[k], G[k] )
        } )
    names( Garg ) <- rep(if (length(t)>0){"-eg"}else{"-g"},length(Garg))
    return(Garg)
}

msarg_M <- function (M,t=numeric(0)) {
    # M is a sparse matrix,
    # so only output the relevant -em commands
    rind <- rowinds(M)
    cind <- colinds(M)
    Marg <- lapply( 1:length(M@x), function (k) {
            c( t, rind[k], cind[k], M@x[k] )
        } )
    names( Marg ) <- rep(if (length(t)>0){"-em"}else{"-m"},length(Marg))
    return(Marg)
}


# The state of a demographic model at a given point in time:
setClass("popArray", representation(
        npop="integer",  # vector of dimensions whose product is number of populations
        N="numeric",     # npop-vector of population sizes
        G="numeric",     # npop-vector of growth rates
        M="Matrix"       # (npop x npop) matrix of migration rates
    ) )

setClass("gridArray", contains="popArray")
setMethod("dim", signature=c(x="gridArray"), definition=function (x) { x@npop } )
setGeneric("layer_inds", function(x,...) { standardGeneric("layer_inds") })
setMethod("layer_inds", signature=c(x="gridArray"), definition=function (x,layer) {
        nrow <- x@npop[1]
        ncol <- x@npop[2]
        return( (layer-1) * nrow * ncol + (1:(nrow*ncol)) )
    } )
setMethod("plot", signature=c(x="gridArray"), definition=function(x,...) {
        # produces nlayers plots
        nlayers <- dim(x)[3]
        layout(matrix(1:(nlayers^2),nrow=nlayers))
        par(oma=c(0.5,0.5,0.5,0.5),mar=c(0.1,0.1,0.1,0.1))
        for (k in 1:nlayers) {
            for (j in 1:nlayers) {
                if (k==j) {
                    plot_layer(x,j,...)
                } else {
                    plot_admixture_layer(x,j,k,...)
                }
            }
        }
    } )

# methods to get row and column indices of nonzero entries of a (sparse) dgCMatrix
setGeneric("rowinds", function(x,...) { standardGeneric("rowinds") } )
setGeneric("colinds", function(x,...) { standardGeneric("colinds") } )
setMethod("rowinds",signature=c(x="dgCMatrix"), definition=function(x) { x@i+1L } )
setMethod("colinds",signature=c(x="dgCMatrix"), definition=function(x) { rep(1:ncol(x),times=diff(x@p)) } )
setMethod("rowinds",signature=c(x="dgTMatrix"), definition=function(x) { x@i+1L } )
setMethod("colinds",signature=c(x="dgTMatrix"), definition=function(x) { x@j+1L } )

plot_admixture_layer <- function (ga,j,k,admix.fac=1,...) {
    # plot patterns of admixture between layers j and k
    # ... admittedly janky.
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    nlayers <- dim(ga)[3]
    rowdummy <- row( Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix(0,nrow=nrow,ncol=ncol) )
    li.j <- layer_inds(ga,j)
    li.k <- layer_inds(ga,k)
    admixmat <- diag( ga@M[li.j,li.k] )
    dim(admixmat) <- c(nrow,ncol)
    plot( row(admixmat), col(admixmat), cex=admix.fac*admixmat, xaxt='n', yaxt='n', xlab='', ylab='', asp=1 )
    return(invisible(admixmat))
}

plot_layer <- function (ga,layer,eps=0.05,lwd.fac=1,length=0.03,...) {
    # draw a picture of a migration matrix
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    rowdummy <- row( Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix(0,nrow=nrow,ncol=ncol) )
    this.M <- ga@M[layers,layers]
    x0 <- rowdummy[rowinds(this.M)]
    x1 <- rowdummy[colinds(this.M)]
    y0 <- coldummy[rowinds(this.M)]
    y1 <- coldummy[colinds(this.M)]
    upshift <- eps*(y1-y0)
    rightshift <- eps*(x1-x0)
    plot(as.vector(rowdummy),as.vector(coldummy),pch=20,xlim=c(0,nrow+1),ylim=c(0,ncol+1),asp=1,xaxt='n',yaxt='n',xlab='',ylab='')
    arrows( x0=x0+upshift, x1=x1+upshift,
       y0=y0+rightshift, y1=y1+rightshift,
       lwd=lwd.fac*this.M@x, length=length, ... )
}


grid_array <- function (nlayers, nrow, ncol, N, mig.rate, admix.rate, G=0) {
    # make a gridArray object
    # consisting of a stack of (nrow x ncol) grids
    #   each having migration rate 'mig.rate' between adjacents
    #   and population sizes N
    # and with correponding pops connected by migration rate admix
    #
    # N and G will be recycled appropriately (in row x col x layer order)
    # mig.rate can be a vector of length nlayers or a single number
    # admix.rate can be a constant or a ( nlayers x nlayers ) matrix
    nn <- nrow*ncol
    N <- rep_len(N,nn*nlayers)
    mig.rate <- rep_len(mig.rate,nlayers)
    admix.rate <- rep_len(admix.rate,nlayers^2)
    dim(admix.rate) <- c(nlayers,nlayers)
    G <- rep_len(G,nn*nlayers)
    # start to build the migration matrix
    M <- kronecker( admix.rate, Diagonal(n=nn,x=1) )
    adj <- grid.adjacency(nrow,ncol,diag=FALSE)
    for (k in 1:nlayers) {
        M[ (k-1)*nn+(1:nn), (k-1)*nn+(1:nn) ] <- mig.rate[k] * adj
    }
    npop <- length(N)
    stopifnot( length(N) == npop && length(G) == npop && ncol(M) == npop && nrow(M) == npop )
    return( new("gridArray",
            npop=as.integer(c(nrow,ncol,nlayers)),
            N=N,
            G=G,
            M=as(M,"dgCMatrix")  # avoids bugs with dgTMatrix subsetting
        ) )
}

modify_grid_layer <- function (ga, layer, dN=1, dG=1, dM=1) {
    # modify a given layer of a grid_array
    # by multiplying N by dN, etcetera.
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    layers <- (layer-1) * nrow * ncol + (1:(nrow*ncol))
    dN <- rep_len(dN,nrow*ncol)
    dG <- rep_len(dG,nrow*ncol)
    dM <- rep_len(dM,(nrow*ncol)^2)
    ga@N[layers] <- ga@N[layers] * dN
    ga@G[layers] <- ga@G[layers] * dG
    ga@M[layers,layers] <- ga@M[layers,layers] * dM
    return(ga)
}

add_barrier <- function (ga, layer, rows=numeric(0), cols=numeric(0)) {
    # remove migration between rows and rows+1 and cols and cols+1
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    dummy <- Matrix(0,nrow=nrow,ncol=ncol)
    for (rn in rows) {
        zeros <- ( ( rowinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(row(dummy)==rn) ) ) & ( colinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(row(dummy)==rn+1) ) ) )
        zeros <- zeros | ( ( rowinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(row(dummy)==rn+1) ) ) & ( colinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(row(dummy)==rn) ) ) )
        ga@M@x[zeros] <- 0
    }
    for (cn in cols) {
        zeros <- ( ( rowinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(col(dummy)==cn) ) ) & ( colinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(col(dummy)==cn+1) ) ) )
        zeros <- zeros | ( ( rowinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(col(dummy)==cn+1) ) ) & ( colinds(ga@M) %in% ( (layer-1)*nrow*ncol + which(col(dummy)==cn) ) ) )
        ga@M@x[zeros] <- 0
    }
    return(ga)
}

modify_migration <- function (ga, layer, doutM=1, dinM=1) {
    # modify inmigration rates on a given layer:
    #   here, dinM is a vector of multiplicative changes to inmigration rates 
    #   here, doutM is a vector of multiplicative changes to outmigration rates 
    # recall that "m ij is the fraction of subpopulation i made up of migrants each generation from subpopulation j."
    #   so that 'inmigration', in proper, forwards, time is m[i,]
    #   and 'outmigration', in proper, forwards, time is m[,j].
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    dinM <- rep_len(dinM,nrow*ncol)
    doutM <- rep_len(doutM,nrow*ncol)
    dM <- outer(dinM,doutM,"*")
    ga@M[layers,layers] <- ga@M[layers,layers] * dM
    return(ga)
}

restrict_patch <- function (ga, layer, xy, r) {
    # zero out mutation rates to/from all populations 
    #   more than distance r from (row,col) location xy
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    rowdummy <- row( Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix(0,nrow=nrow,ncol=ncol) )
    base.ind <- which( ( rowdummy == xy[1] ) & ( coldummy == xy[2] ) )
    gdists <- sqrt( (rowdummy-xy[1])^2 + (coldummy-xy[2])^2 )
    zero.these <- layers[ gdists > r ]
    ga@M[zero.these,] <- 0.0
    ga@M[,zero.these] <- 0.0
    return(ga)
}


##
# stuff to make a grid matrix

grid.adjacency <- function (nrow,ncol=nrow,diag=TRUE,symmetric=TRUE) {
    # mappings between index in a matrix of height (in rows) nrow
    #   ZERO-BASED: (i,j) and column-oriented (k)
    #   ALLOW indices outside the grid
    .ob <- function (ij,nrow,ncol){ ( ij[,1] >= 0 ) & ( ij[,1] < nrow ) & ( ij[,2] >= 0 ) & ( ij[,2] < ncol ) } # outside grid
    ij.to.k <- function (ij,nrow,ncol) { if (is.null(dim(ij))) { dim(ij) <- c(1,length(ij)) }; ifelse( .ob(ij,nrow,ncol), ij[,1,drop=FALSE]+ij[,2,drop=FALSE]*nrow, NA ) }
    k.to.ij <- function (k,nrow) { cbind( k%%nrow, k%/%nrow ) }
    shift <- function (dij,k,nrow,ncol) { ij.to.k( sweep( k.to.ij(k,nrow), 2, as.integer(dij), "+" ), nrow, ncol ) }
    # returns adjacency matrix for a grid of height nrow and width ncol
    nn <- 0:(nrow*ncol-1)
    nreps <- if(diag){5}else{4}
    adj <- data.frame(
            i=rep(nn, nreps),
            j=c( if(diag){nn}else{NULL},
                 shift(c(+1,0),nn,nrow,ncol),
                 shift(c(-1,0),nn,nrow,ncol),
                 shift(c(0,+1),nn,nrow,ncol),
                 shift(c(0,-1),nn,nrow,ncol)
                 )
            )
    # on the boundary?
    usethese <- ! ( is.na(adj$j) | !.ob(k.to.ij(adj$j,nrow),nrow,ncol) | !.ob(k.to.ij(adj$i,nrow),nrow,ncol) )
    if (symmetric) { usethese <- ( usethese & ( adj$i <= adj$j ) ) }
    # add 1 since we worked 0-based above; sparseMatrix (unlike underlying representation) is 1-based.
    A <- with( subset(adj, usethese ), sparseMatrix( i=i+1L, j=j+1L, x=1.0, dims=c(nrow*ncol,nrow*ncol), symmetric=symmetric ) )
    return(A)
}


