
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

setMethod("dim", signature=c(x="demography"), definition=function (x) { dim(x[[1]]) } )
setMethod("length", signature=c(x="demography"), definition=function (x) { length(x@t) } )
setMethod("[[", signature=c(x="demography",i="ANY"), definition=function (x,i) { x@popStates[[i]] } )
setMethod("[[<-", signature=c(x="demography",i="ANY",value="ANY"), definition=function (x,i,value) { x@popStates[[i]]<-value; return(x) } )

add_to_demography <- function (dem,dt,fn=identity,pop,...,tnew=dt+if(length(dem@t)>0){dem@t[length(dem@t)]}else{0}) {
    # add to a demography by applying a modifier function fn 
    #  to another population model,
    #  by default the previous state
    if (inherits(dem,"popArray")) {
        dem <- new("demography",popStates=list(dem),t=numeric(0))
    }
    nt <- length(dem@popStates)
    if (missing(pop)) { pop <- dem@popStates[[nt]] }
    dem@popStates <- c( dem@popStates, fn(pop,...) )
    stopifnot( tnew >= max(c(0,dem@t)) )
    dem@t <- c(dem@t,tnew)
    return(dem)
}

demography <- function (ga,t) { new("demography",popStates=list(ga),t=numeric(0)) }

### Stuff for writing ms arguments
setClass("msarg", contains="namedList")
setMethod("toString", signature=c(x="msarg"), definition=function (x,sep="\n") {
        text.x <-  paste( paste( paste(names(x),sapply(lapply(x,unlist),paste,collapse=' ')), collapse=sep ), "\n", sep='' )
        return(invisible(text.x))
    } )
setMethod("print", signature=c(x="msarg"), definition=function (x,sep="\n",...) {
        text.x <-  toString(x,sep=sep)
        cat( text.x, ... )
        return(invisible(text.x))
    } )

setGeneric("msarg", function(x,...) { standardGeneric("msarg") })
setMethod("msarg", signature=c(x="demography"), definition=function (x,nsamp,scale.migration=TRUE,...) {
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
        if (NCOL(nsamp)==4) {
            # assume this is in sample.config format
            nsamp <- nsample_vector(x[[1]],nsamp)
        }
        stopifnot( length(nsamp) ==  prod(dim(x[[1]])) )
        arglist <- c( list(
                        "-I"=c( prod(dim(x[[1]])), nsamp )
                        ),
                    msarg_N(x[[1]]),  # note that N must come before G, as N sets G to zero.
                    msarg_G(x[[1]]),
                    msarg_M(x[[1]],scale.migration=scale.migration)
                )
        for (k in seq_along(x@t)) {
            t <- x@t[k]
            arglist <- c( 
                    arglist,
                    msarg_N(x[[k+1]],t,previous=x[[k]]),  # note that N must come before G, as N sets G to zero.
                    msarg_G(x[[k+1]],t,previous=x[[k]]),
                    msarg_M(x[[k+1]],t,previous=x[[k]],scale.migration=scale.migration)
                )
        }
        return(new("msarg",.Data=arglist,names=names(arglist)))
    } )

msarg_N <- function (ga,t=numeric(0),previous,eps=min(ga@N[ga@N>0])/1000) {
    # Ns that haven't changed don't matter
    #   and replace zeros with eps
    N <- ga@N
    N[N<=0] <- eps
    # so only output the relevant -eg commands
    dothese <- if (missing(previous)) { seq_along(N) } else { which(N!=previous@N) }
    Narg <- lapply( seq_along(dothese), function (k) {
            c( t, dothese[k], N[k] )
        } )
    names( Narg ) <- rep(if (length(t)>0){"-en"}else{"-n"},length(Narg))
    return(Narg)
}

msarg_G <- function (ga,t=numeric(0),previous) {
    # Gs that haven't changed don't matter
    #   so only output the relevant -eg commands
    # This could be more efficient, omitting zeros at the start, etc.
    G <- ga@G
    dothese <- if (missing(previous)) { seq_along(G) } else { which(G!=previous@G) }
    Garg <- lapply( seq_along(dothese), function (k) {
            c( t, dothese[k], G[k] )
        } )
    names( Garg ) <- rep(if (length(t)>0){"-eg"}else{"-g"},length(Garg))
    return(Garg)
}

msarg_M <- function (ga,t=numeric(0),previous,scale.migration=TRUE) {
    # M is a sparse matrix,
    #   so only output the relevant -em commands
    # also: need to transpose and scale by population sizes
    #   since M is forwards-times migration rates,
    #   and ms takes revese-time ("lineage") migration rates.
    # (optionally skip this step to allow ms-style migration input)
    M <- if (scale.migration) { t(lineage_M(ga)) } else { ga@M }
    rind <- rowinds(M)  # yes this is inefficient
    cind <- colinds(M)
    dothese <- if (missing(previous)) {
            seq_along(M@x)
        } else {
            old.M <- lineage_M(previous)
            which( (M@x != old.M@x) | (rind!=rowinds(old.M)) | (cind!=colinds(old.M)) )
        }
    Marg <- lapply( dothese, function (k) {
            c( t, rind[k], cind[k], M@x[k] )
        } )
    names( Marg ) <- rep(if (length(t)>0){"-em"}else{"-m"},length(Marg))
    return(Marg)
}

lineage_M <- function (ga,eps=min(ga@N[ga@N>0])/1000) {
    M <- ga@M
    rind <- rowinds(M)
    cind <- colinds(M)
    M@x <- M@x * (ga@N[rind]+eps) / (ga@N[cind]+eps)  # scale by pop sizes
    return(M)
}

### stuff for calling ms
# general strategy is: put everything in a subdirectory

setGeneric("run_ms", function(x,...) { standardGeneric("run_ms") })
setMethod("run_ms", signature=c(x="popArray"), definition=function (x,...) { run_ms(demography(x),...) } )
setMethod("run_ms", signature=c(x="demography"), definition=function (x,nsamp,outdir,theta,trees=FALSE,nreps=1,tofile=TRUE,...) {
        if (NCOL(nsamp)==4) {
            # assume this is in sample.config format
            nsamp <- nsample_vector(x[[1]],nsamp)
        }
        msarg <- msarg(x,nsamp,...)
        if (tofile) {
            if (missing(outdir)) {
                while (TRUE) {
                    outdir <- paste( "msarg", formatC( 1e6*runif(1), width=6, format='d', flag='0' ), sep="_" )
                    if (!file.exists(outdir)) { break }
                }
            }
            cat("Outputting to ", outdir, " .\n")
            dir.create(outdir,showWarnings=FALSE)
            msarg.file <- file.path(outdir,"msarg.txt")
            print(msarg,file=msarg.file,sep='  ')
        }
        ms.call <- paste( "ms", sum(nsamp), nreps, 
            if (!missing(theta)) { paste("-t",theta) } else {""},
            if (trees) { "-T" } else {""},
            if (tofile) { paste("-f", msarg.file) } else { toString(msarg,sep=' ') },
            if (tofile) { paste( ">", file.path(outdir,"msoutput.txt") ) } else { "" } 
            ) 
        ms.results <- system( ms.call, intern=TRUE )
        return( if (tofile) { invisible(outdir) } else { ms.results } )
    } )

trees_from_ms <- function (ms.output) {
    # extract the trees from ms.output
    if ( (length(ms.output)==1) && (file.exists(ms.output)) ) {
        ms.output <- scan(ms.output,what='char')
    }
    almost <- which( grepl("^//",ms.output) )
    lapply( ms.output[almost+1], function (x) { ape::read.tree(text=x) } )
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
setGeneric("nlayers", function(x,...) { standardGeneric("nlayers") })
setMethod("nlayers", signature=c(x="gridArray"), definition=function (x) { dim(x)[3] } )
setGeneric("layer_inds", function(x,...) { standardGeneric("layer_inds") })
setMethod("layer_inds", signature=c(x="gridArray"), definition=function (x,layer) {
        nrow <- x@npop[1]
        ncol <- x@npop[2]
        return( (layer-1) * nrow * ncol + (1:(nrow*ncol)) )
    } )
setMethod("plot", signature=c(x="gridArray"), definition=function(x,layers=seq_len(nlayers(x)),do.layout=TRUE,...) {
        # produces nlayers plots
        nlayers <- nlayers(x)
        if (do.layout) {
            layout(matrix(1:(nlayers^2),nrow=nlayers))
            opar <- par(oma=c(0.5,0.5,0.5,0.5),mar=c(0.1,0.1,0.1,0.1))
        }
        for (k in layers) {
            for (j in layers) {
                if (k==j) {
                    plot_layer(x,j,...)
                } else {
                    plot_admixture_layer(x,j,k,...)
                }
            }
        }
        if (do.layout) { par(opar) }
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
    admixmat <- diag( x=ga@M[li.j,li.k], nrow=length(li.j) )
    dim(admixmat) <- c(nrow,ncol)
    plot( row(admixmat), col(admixmat), cex=admix.fac*admixmat, xaxt='n', yaxt='n', xlab='', ylab='', asp=1 )
    return(invisible(admixmat))
}

plot_layer <- function (ga,layer,eps=0.05,lwd.fac=1,cex.fac=2/quantile(ga@N,.9),length=0.03,alpha=0.05,N.eps=1e-3,...) {
    # draw a picture of a migration matrix
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    lineage.M <- lineage_M(ga)
    # for plotting purposes, remove lineage migration rates between 'zeroed' populations.
    if (any(ga@N<=N.eps)) {
        lineage.M[ga@N<=N.eps,] <- 0.0
        lineage.M[,ga@N<=N.eps] <- 0.0  
    }
    if (nrow==1 && ncol==1) {
        plot( 0, type='n', xlab='', ylab='', xaxt='n', yaxt='n' )
    } else {
        rowdummy <- row( Matrix(0,nrow=nrow,ncol=ncol) )
        coldummy <- col( Matrix(0,nrow=nrow,ncol=ncol) )
        this.M <- ga@M[layers,layers]
        rind <- rowinds(this.M)
        cind <- colinds(this.M)
        this.lin.M <- lineage.M[layers,layers][cbind(rind,cind)]
        x0 <- rowdummy[rind]
        x1 <- rowdummy[cind]
        y0 <- coldummy[rind]
        y1 <- coldummy[cind]
        upshift <- eps*(y1-y0)
        rightshift <- eps*(x1-x0)
        plot(
            as.vector(rowdummy),as.vector(coldummy),
            cex=cex.fac*ga@N[layers],
            pch=20,xlim=c(0,nrow+1),ylim=c(0,ncol+1),
            asp=1,xaxt='n',yaxt='n',xlab='',ylab='')
        arrows( x0=x0+upshift, x1=x1+upshift,
            y0=y0+rightshift, y1=y1+rightshift,
            col=sapply((this.lin.M/max(lineage.M@x,na.rm=TRUE))^alpha,function(a){adjustcolor("black",a)}),
            lwd=lwd.fac*this.M@x, length=length, ... )
    }
}


grid_array <- function (nlayers=1, nrow, ncol, N=1, mig.rate=0, admix.rate=0, G=0) {
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
    # avoid bugs with dgTMatrix subsetting
    if (class(M)=="dtTMatrix") { # no method for going from dtTMatrix to dgCMatrix??
        M <- as(M,"dgTMatrix")
    }
    M <- as(M,"dgCMatrix")
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
            M=M
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
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    dinM <- rep_len(dinM,nrow*ncol)
    doutM <- rep_len(doutM,nrow*ncol)
    dM <- outer(doutM,dinM,"*")
    ga@M[layers,layers] <- ga@M[layers,layers] * dM
    return(ga)
}

restrict_patch <- function (ga, layer, xy, r) {
    # zero out all populations
    #   more than distance r from (row,col) location xy
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    rowdummy <- row( Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix(0,nrow=nrow,ncol=ncol) )
    base.ind <- which( ( rowdummy == xy[1] ) & ( coldummy == xy[2] ) )
    gdists <- sqrt( (rowdummy-xy[1])^2 + (coldummy-xy[2])^2 )
    zero.these <- layers[ gdists > r ]
    ga@N[zero.these] <- 0
    return(ga)
}

##
# stuff for selecting sample locations
# format:
#  row column layer n
#   1    1      1   10
#   2    4      1   12
#   ....

plot_sample_config <- function ( dem, sample.config, sample.cols=rainbow(nrow(sample.config)), add=FALSE, ... ) {
    # first put in the same order as processed by ms:
    sample.config <- sample.config[do.call(order,list(sample.config[,1],sample.config[,2],sample.config[,3])),]
    if (!add) { 
        plot( sample.config[,1:2], type='n', xlim=c(0,nrow(dem)+1), ylim=c(0,ncol(dem)+1), xlab='', ylab='', xaxt='n', yaxt='n', asp=1, bty='n' ) 
        rect( xleft=1, ybottom=1, xright=nrow(dem), ytop=ncol(dem) )
    }
    text( sample.config[,1:2], 
            labels=tapply( seq_len(sum(sample.config[,4])), rep(1:nrow(sample.config),sample.config[,4]), paste, collapse=',' ),
            col=sample.cols )
}

nsample_vector <- function (ga, samps, dims=dim(ga)) {
    # samps is a four-column matrix,
    #    of row number, column number, layer number, and number of samples.
    # returns a vector with length equal to the number of demes in ga
    #    that gives the number of samples in each  (this is probably quite long!)
    if (nrow(samps) != nrow(unique(samps))) {
        stop("Sampling locations are not unique.")
    }
    sample.config <- array(0,dim=dims)
    sample.config[ as.matrix(samps)[,1:3] ] <- samps[,4]
    return( as.vector(sample.config) )
}

sample_locations <- function (ga, n, each=1, dims=dim(ga)) {
    # Sample uniformly at random, with replacement, `n` locations to sample from.
    n <- each * tabulate( sample.int( prod(dims), n, replace=TRUE ), nbins=prod(dims) )
    x <- arrayInd( which(n>0), .dim=dims )
    colnames(x) <- c("row","col","layer")
    return( cbind(x,n=n[n>0]) )
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


