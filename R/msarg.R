

##
# A demographic model consists of:
#  - a set of populations (demes)
#  - their sizes
#  - and their pairwise migration rates
# all possibly varying in time, and so hence also:
#  - an identification of populations between time points

#' @importClassesFrom Matrix dgCMatrix dgTMatrix Matrix
NULL

#' @importClassesFrom Biostrings DNAStringSet
NULL

# An entire demographic model strings together states:

#' @export
methods::setClass("ms_demog", representation(
        popStates="list",       # a list of popArrays
        t="numeric"             # vector of lengths of times ago of transitions between the states
    ) )

# The state of a demographic model at a given point in time:

#' @export
methods::setClass("popArray", representation(
        npop="integer",  # vector of dimensions whose product is number of populations
        N="numeric",     # npop-vector of population sizes
        G="numeric",     # npop-vector of growth rates
        M="Matrix"       # (npop x npop) matrix of migration rates
    ) )

#' @export
methods::setClass("gridArray", contains="popArray")

methods::setMethod("dim", signature=c(x="ms_demog"), definition=function (x) { dim(x[[1]]) } )
methods::setMethod("length", signature=c(x="ms_demog"), definition=function (x) { length(x@popStates) } )

#' @export
methods::setAs("ms_demog", "list", def=function (from) { from@popStates } )
methods::setMethod("[[", signature=c(x="ms_demog",i="ANY"), definition=function (x,i) { x@popStates[[i]] } )
methods::setMethod("[[<-", signature=c(x="ms_demog",i="ANY",value="ANY"), definition=function (x,i,value) { x@popStates[[i]]<-value; return(x) } )
methods::setMethod("plot", signature=c(x="ms_demog"), definition=function (x,...) {
        NN <- sapply( x@popStates, slot, "N" )
        cex.fac <- 5/max(1,NN)
        for (k in seq_along(x@popStates)) {
            readline(paste("Hit enter for t = ", c(0,x@t)[k]))
            plot(x[[k]], cex.fac=cex.fac, ...)
        }
    } )

#' @export
ms_demog <- function (ga,t) { new("ms_demog",popStates=list(ga),t=numeric(0)) }

### Stuff for writing ms arguments

#' @export
methods::setClass("msarg", contains="namedList")
methods::setMethod("toString", signature=c(x="msarg"), definition=function (x,sep="\n") {
        text.x <-  paste( paste( paste(names(x),sapply(lapply(x,unlist),paste,collapse=' ')), collapse=sep ), "\n", sep='' )
        return(invisible(text.x))
    } )
methods::setMethod("print", signature=c(x="msarg"), definition=function (x,sep="\n",...) {
        text.x <-  toString(x,sep=sep)
        cat( text.x, ... )
        return(invisible(text.x))
    } )

#' Construct an ms Command Line From a Demography
#'
#' @param x The demography object.
#' @param nsamp The sample configuration.
#' @param scale.migration Given migration rates are forwards-time (so must rescale to reverse-time rates for ms).
#' @param N.eps Set populations with zero population size to this value, to avoid noncommunicating subpopulations.
#' @param ... Other parameters passed to \code{ms}.
#' @export msarg
#' @name msarg
methods::setGeneric("msarg", function(x,...) { standardGeneric("msarg") })

methods::setMethod("msarg", signature=c(x="ms_demog"), 
   definition=function (x,
                        nsamp,
                        scale.migration=TRUE,
                        N.eps=NULL,
                        ...) {
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
        if (prod(dim(x[[1]]))==1) {
            # ms won't allow -I with only one pop
            stop("Sorry, can't do single populations at the moment.")
        }
        more.args <- list(...)
        if (length(more.args)>0) { names(more.args) <- paste("-",names(more.args),sep='') }
        arglist <- c( 
                    more.args,  # can pass in other things, like -p n or -seeds x1 x2 x3
                    list(
                        "-I"=c( prod(dim(x[[1]])), nsamp )
                        ),
                    msarg_N(x[[1]],eps=N.eps),  # note that N must come before G, as N sets G to zero.
                    msarg_G(x[[1]]),
                    msarg_M(x[[1]],scale.migration=scale.migration,eps=N.eps)
                )
        for (k in seq_along(x@t)) {
            t <- x@t[k]
            arglist <- c( 
                    arglist,
                    msarg_N(x[[k+1]],t,previous=x[[k]],eps=N.eps),  # note that N must come before G, as N sets G to zero.
                    msarg_G(x[[k+1]],t,previous=x[[k]]),
                    msarg_M(x[[k+1]],t,previous=x[[k]],scale.migration=scale.migration,eps=N.eps)
                )
        }
        return(new("msarg",.Data=arglist,names=names(arglist)))
    } )

msarg_N <- function (ga,t=numeric(0),previous,eps) {
    # Ns that haven't changed don't matter
    #   and replace zeros with eps
    if (all(ga@N==0)) { stop(sprintf("All population sizes are less than %f.",eps)) }
    if (missing(eps) || is.null(eps)) { eps <- min(ga@N[ga@N>0])/1000 }
    N <- ga@N
    N[N<eps] <- eps
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

msarg_M <- function (ga,t=numeric(0),previous,scale.migration=TRUE,eps=NULL) {
    # M is a sparse matrix,
    #   so only output the relevant -em commands
    # also: need to transpose and scale by population sizes
    #   since M is forwards-times migration rates,
    #   and ms takes reverse-time ("lineage") migration rates.
    # (optionally skip this step to allow ms-style migration input)
    M <- if (scale.migration) { Matrix::t(lineage_M(ga,eps=eps)) } else { ga@M }
    rind <- rowinds(M)  # yes this is inefficient
    cind <- colinds(M)
    # don't do the diagonal
    dothese <- ( rind != cind )
    if (any(M@x[dothese]<0)) { stop("M has negative entries.") }
    # don't respecify migration rates that haven't changed
    if (!missing(previous)) {
        # note: if structure of matrix has changed will respecify some that don't need to be (no harm)
        old.M <- lineage_M(previous)
        if ( (length(M@x)==length(old.M@x)) ) {
            dothese <- ( dothese & ( (M@x != old.M@x) | (rind!=rowinds(old.M)) | (cind!=colinds(old.M)) ) )
        }
    }
    Marg <- lapply( which(dothese), function (k) {
            c( t, rind[k], cind[k], M@x[k] )
        } )
    names( Marg ) <- rep(if (length(t)>0){"-em"}else{"-m"},length(Marg))
    return(Marg)
}

#' @export
lineage_M <- function (ga,eps) {
    if (missing(eps) || is.null(eps)) { eps <- min(ga@N[ga@N>0])/1000 }
    M <- ga@M
    rind <- rowinds(M)
    cind <- colinds(M)
    M@x <- M@x * (ga@N[rind]+eps) / (ga@N[cind]+eps)  # scale by pop sizes
    return(M)
}

### stuff for calling ms
# general strategy is: put everything in a subdirectory

#' @export run_ms
methods::setGeneric("run_ms", function(x,...) { standardGeneric("run_ms") })
methods::setMethod("run_ms", signature=c(x="popArray"), definition=function (x,...) { run_ms(ms_demog(x),...) } )
methods::setMethod("run_ms", signature=c(x="ms_demog"), 
    definition=function (x,
                         nsamp,
                         outdir,
                         theta,
                         trees=FALSE,
                         nreps=1,
                         tofile=TRUE,
                         ms.binary="ms",
                         ...) {
        if (missing(theta) & !trees) { stop("Must specify either theta or trees=TRUE.") }
        if (NCOL(nsamp)==4) {
            # assume this is in sample.config format
            nsamp <- nsample_vector(x[[1]],nsamp)
        }
        ms.arg <- msarg(x,nsamp,...)
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
            cat(toString(ms.arg),file=msarg.file)
        }
        ms.call <- paste( ms.binary, sum(nsamp), nreps, 
            if (!missing(theta)) { paste("-t",theta) } else {""},
            if (trees) { "-T" } else {""},
            if (tofile) { paste("-f", msarg.file) } else { toString(ms.arg,sep=' ') },
            if (tofile) { paste( ">", file.path(outdir,"msoutput.txt") ) } else { "" } 
            ) 
        ms.results <- system( ms.call, intern=TRUE )
        return( if (tofile) { invisible(outdir) } else { ms.results } )
    } )

##
# processing tree output

#' @export
trees_from_ms <- function (ms.output) {
    # extract the trees from ms.output
    if ( (length(ms.output)==1) && (file.exists(ms.output)) ) {
        ms.output <- scan(ms.output,what='char')
    }
    almost <- which( grepl("^//",ms.output) )
    out <- lapply( ms.output[almost+1], function (x) { 
               tryCatch( ape::read.tree(text=x),
                        error=function(e) { 
                            # cat("Warning: bad tree.\n")
                            # cat("  ", e$message, "\n")
                            NULL
                        } )
        } )
    skip <- sapply( out, is.null )
    if (sum(skip)>0) { warning(sprintf("skipping %d bad trees.",sum(skip))) }
    return( out[!skip] )
}

#' @export
tree_dists <- function (trees,sample.config) {
    # get matrix of tree distances
    #  **in the right order**
    tip_order <- tip_order_fn(sample.config)
    lapply( trees, function (tree) {
            ord <- tip_order(tree)
            ape::cophenetic.phylo(tree)[ord,ord]
        } )
}

#' @export
tip_order_fn <- function (sample.config) {
    # return the vector of indices that will reorder the tips of a tree to match sample.config
    # first put in the same order as processed by ms:
    sample.config.ord <- order(sample.config[,3],sample.config[,2],sample.config[,1])
    #  after reordering, sample.config[k,] is the k-th *group* of output*s* from ms
    sample.config <- sample.config[sample.config.ord,,drop=FALSE]
    return( function (tree) {
        # tip.ord[k] says which *reordered* group the k-th tip is in
        tip.ord <- findInterval(as.numeric(tree$tip.label),1+c(0,cumsum(sample.config[,4])))
        # so the k-th tip is in the sample.config.ord[tip.ord[k]]-th original group
        # rank(sample.config.ord[tip.ord],ties.method="first")  # this MISTAKE ends up putting nearby leaves close to each other, strangely
        order(sample.config.ord[tip.ord],as.numeric(tree$tip.label))
    } )
}


# draw vertical lines at particular times on a phylogeny
#' @export
abline_phylo <- function (v=NULL,
                          backward=TRUE,
                          root.time=NULL,
                          lastPP= get("last_plot.phylo", envir = .PlotPhyloEnv),
                          ...) {
    # from ape::axisPhylo
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards")) range(lastPP$xx) else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp)) tscale <- tscale[2:1]
    if (!is.null(root.time)) {
        tscale <- tscale + root.time
        if (backward) tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    abline(v=beta*v+alpha,...)
    return(invisible(c(alpha=alpha,beta=beta)))
}

## methods for gridArray's
methods::setMethod("dim", signature=c(x="gridArray"), definition=function (x) { x@npop } )

#' @export nlayers
methods::setGeneric("nlayers", function(x,...) { standardGeneric("nlayers") })
methods::setMethod("nlayers", signature=c(x="gridArray"), definition=function (x) { dim(x)[3] } )

#' @export layer_inds
methods::setGeneric("layer_inds", function(x,...) { standardGeneric("layer_inds") })
methods::setMethod("layer_inds", signature=c(x="gridArray"), definition=function (x,layer) {
        # tells you which indices in e.g. x@N correspond to layer number 'layer'
        nrow <- x@npop[1]
        ncol <- x@npop[2]
        return( (layer-1) * nrow * ncol + (1:(nrow*ncol)) )
    } )
methods::setMethod("plot", signature=c(x="gridArray"), definition=function(x,layers=seq_len(nlayers(x)),do.layout=TRUE,...) {
        # produces nlayers plots
        nlayers <- nlayers(x)
        if (do.layout) {
            layout(matrix(1:(nlayers^2),nrow=nlayers))
            opar <- par(oma=c(0.5,0.5,0.5,0.5),mar=c(0.1,0.1,0.1,0.1))
            on.exit( par(opar), add=TRUE )
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
    } )

# methods to get row and column indices of nonzero entries of a (sparse) dgCMatrix

#' @export rowinds
methods::setGeneric("rowinds", function(x,...) { standardGeneric("rowinds") } )

#' @export colinds
methods::setGeneric("colinds", function(x,...) { standardGeneric("colinds") } )
methods::setMethod("rowinds",signature=c(x="dgCMatrix"), definition=function(x) { x@i+1L } )
methods::setMethod("colinds",signature=c(x="dgCMatrix"), definition=function(x) { rep(1:ncol(x),times=diff(x@p)) } )
methods::setMethod("rowinds",signature=c(x="dgTMatrix"), definition=function(x) { x@i+1L } )
methods::setMethod("colinds",signature=c(x="dgTMatrix"), definition=function(x) { x@j+1L } )

#' @export
plot_admixture_layer <- function (ga,j,k,admix.fac=1,...) {
    # plot patterns of admixture between layers j and k
    # ... admittedly janky.
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    nlayers <- dim(ga)[3]
    rowdummy <- row( Matrix::Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix::Matrix(0,nrow=nrow,ncol=ncol) )
    li.j <- layer_inds(ga,j)
    li.k <- layer_inds(ga,k)
    admixmat <- diag( x=ga@M[li.j,li.k], nrow=length(li.j) )
    dim(admixmat) <- c(nrow,ncol)
    plot( row(admixmat), col(admixmat), cex=admix.fac*admixmat, xaxt='n', yaxt='n', xlab='', ylab='', asp=1 )
    return(invisible(admixmat))
}

.exlims <- function (x) { 1.05 * (x-mean(x)) + mean(x) }

#' @export
plot_layer <- function (ga,
                        layer,
                        eps=0.05,
                        lwd.fac=1,
                        cex.fac=2/quantile(pmax(10,ga@N[ga@N>0]),.9),
                        length=0.03,
                        alpha=0.05,
                        N.eps=1e-3,
                        do.arrows=TRUE,
                        arrow.direction=c("forwards","reverse"), # do arrows for forwards-time migration or reverse?
                        add=FALSE,
                        xy=cbind( x=as.vector(row( Matrix::Matrix(0,nrow=nrow(ga),ncol=ncol(ga)) )),
                                  y=as.vector(col( Matrix::Matrix(0,nrow=nrow(ga),ncol=ncol(ga)) )) ),
                          # location of the points on the map
                        ...) {
    if (all(ga@N==0)) { stop("All population sizes are zero.") }
    # draw a picture of a migration matrix
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    if (nrow==1 && ncol==1) {
        if (!add) { plot( 0, type='n', xlab='', ylab='', xaxt='n', yaxt='n' ) }
    } else {
        if (!add) { 
            plot( xy, type='n',
                    xlim=.exlims(range(xy[,"x"],finite=TRUE)),
                    ylim=.exlims(range(xy[,"y"],finite=TRUE)),
                    asp=1, xaxt='n', yaxt='n', xlab='', ylab='')
        }
        if (do.arrows) {
            lineage.M <- lineage_M(ga)
            this.M <- ga@M[layers,layers]
            # for plotting purposes, 
            # remove lineage migration rates between 'zeroed' populations.
            remove.links <- ( is.na(ga@N) | ( ga@N<=N.eps ) )
            if (any(remove.links)) {
                lineage.M[remove.links,] <- 0.0
                lineage.M[,remove.links] <- 0.0
                this.M[remove.links,] <- 0.0
                this.M[,remove.links] <- 0.0
            }
            # indices corresponding to rows and columns of this.M@x
            rind <- rowinds(this.M)
            cind <- colinds(this.M)
            use.these <- ( (rind != cind) & (this.M@x > 0) )
            rind <- rind[use.these]
            cind <- cind[use.these]
            nmags <- 64
            arrow.magnitudes <- cut( switch(match.arg(arrow.direction),
                        "forwards" = this.M[layers,layers][cbind(rind,cind)],
                        "reverse" = lineage.M[layers,layers][cbind(rind,cind)]
                    ), nmags )
            # this.lin.M <- lineage.M[layers,layers][cbind(rind,cind)]
            # this.M.vals <- this.M[layers,layers][cbind(rind,cind)]
            x0 <- xy[rind,"x"]
            x1 <- xy[cind,"x"]
            y0 <- xy[rind,"y"]
            y1 <- xy[cind,"y"]
            upshift <- eps*(y1-y0)
            rightshift <- eps*(x1-x0)
            migr.col <- rev(heat.colors(1.5*nmags)[1:nmags])
            arrows( x0=x0+upshift, x1=x1+upshift,
                    y0=y0+rightshift, y1=y1+rightshift,
                    # col=sapply((this.lin.M/max(lineage.M@x,na.rm=TRUE))^alpha,function(a){adjustcolor("black",a)}),
                    # lwd=lwd.fac*this.lin.M, length=length, ... )
                    col=migr.col[arrow.magnitudes], length=length, ... )
        }
        points( xy, cex=cex.fac*ga@N[layers], pch=20 )
    }
}

###
# stuff for manipulating grids and demographies

#' @export
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
    M <- kronecker( admix.rate, Matrix::Diagonal(n=nn,x=1) )
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

#' @export
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

#' @export
add_barrier <- function (ga, layer, rows=numeric(0), cols=numeric(0)) {
    # remove migration between rows and rows+1 and cols and cols+1
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    dummy <- Matrix::Matrix(0,nrow=nrow,ncol=ncol)
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

#' @export
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

#' @export
restrict_patch <- function (ga, layer, xy, r) {
    # zero out all populations
    #   more than distance r from (row,col) location xy
    layers <- layer_inds(ga,layer)
    nrow <- nrow(ga)
    ncol <- ncol(ga)
    rowdummy <- row( Matrix::Matrix(0,nrow=nrow,ncol=ncol) )
    coldummy <- col( Matrix::Matrix(0,nrow=nrow,ncol=ncol) )
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

## TO-DO: add a class for sample config

#' @export
sort_sample_config <- function (sample.config) {
    sample.config[order(sample.config[,3],sample.config[,2],sample.config[,1]),,drop=FALSE]
}

#' @export
sample_config_index <- function( sample.config, ga, .dim=dim(ga) ) {
    # from the row,col,layer information in sample.config,
    # obtain the index of the sampled locations
    .dim <- c(.dim,rep(1,3-length(.dim)))
    dummy <- array(seq_len(prod(.dim)), dim=.dim)
    return( dummy[ sample.config[,1:3] ] )
}

#' @export
plot_sample_config <- function ( 
                        dem,
                        sample.config,
                        col=rainbow(nrow(sample.config)),
                        add=FALSE,
                        xy=cbind( x=as.vector(row( Matrix::Matrix(0,nrow=nrow(dem),ncol=ncol(dem)) )),
                                  y=as.vector(col( Matrix::Matrix(0,nrow=nrow(dem),ncol=ncol(dem)) )) ),
                        labels=tapply( seq_len(sum(sample.config[,4])), rep(1:nrow(sample.config),sample.config[,4]), paste, collapse=',' ),
                        ... ) {
    # first put in the same order as processed by ms:
    sample.config <- sort_sample_config(sample.config)
    sample.xy <- xy[sample_config_index( sample.config, dem ),]
    if (!add) { 
        # plot( sample.config[,1:2], type='n', xlim=c(0,nrow(dem)+1), ylim=c(0,ncol(dem)+1), xlab='', ylab='', xaxt='n', yaxt='n', asp=1, bty='n' ) 
        plot( sample.xy, type='n', xlim=.exlims(range(xy[,"x"])), ylim=.exlims(range(xy[,"y"])), 
             xlab='', ylab='', xaxt='n', yaxt='n', asp=1, bty='n' ) 
        rect( xleft=min(xy[,"x"]), ybottom=min(xy[,"y"]), xright=max(xy[,"x"]), ytop=max(xy[,"y"]) )
    }
    text( sample.xy, labels=labels, col=col, ... )
}

#' @export
nsample_vector <- function (ga, samps, dims=dim(ga)) {
    # samps is a four-column matrix,
    #    of row number, column number, layer number, and number of samples.
    # returns a vector with length equal to the number of demes in ga
    #    that gives the number of samples in each  (this is probably quite long!)
    if (nrow(samps) != nrow(unique(samps))) {
        stop("Sampling locations are not unique.")
    }
    if ( any( samps[,1] > dims[1] ) | any( samps[,2] > dims[2] ) | any( samps[,3] > dims[3] ) ) {
        stop("Sampling locations are out of bounds.")
    }
    if (any(diff(order(samps[,3],samps[,2],samps[,1]))<0)) {
        warning("sample configuration not in sorted order: output from ms will not be in the same order as the configuration.")
    }
    sample.config <- array(0,dim=dims)
    sample.config[ as.matrix(samps)[,1:3] ] <- samps[,4]
    return( as.vector(sample.config) )
}

#' @export
sample_locations <- function (ga, n, each=1, dims=dim(ga)) {
    # Sample uniformly at random, with replacement, `n` locations to sample from.
    n <- each * tabulate( sample.int( prod(dims), n, replace=TRUE ), nbins=prod(dims) )
    x <- arrayInd( which(n>0), .dim=dims )
    colnames(x) <- c("row","col","layer")
    return( cbind(x,n=n[n>0]) )
}

#' @export
sample_real_locs <- function (real.locs, ga, each=1, dims=dim(ga)) {
    # Translate "real" locations to their nearest grid point
    locs <- cbind( i=1+floor(dim(ga)[1]*real.locs[,1]), j=1+floor(dim(ga)[2]*real.locs[,2]) )
    n <- each * tabulate( (locs[,"i"]-1)*dim(ga)[1]+locs[,"j"], nbins=prod(dims) )
    x <- arrayInd( which(n>0), .dim=dims )
    colnames(x) <- c("row","col","layer")
    sample.config <- cbind(x,n=n[n>0]) 
    # output in the order returned by ms
    sample.config <- sort_sample_config(sample.config)
    return(sample.config)
}

#' @export
distance_from_sample <- function (sample.config,rowscale=1,colscale=1) {
    # return a matrix of distances between all samples
    x <- sample.config[,1]*rowscale
    y <- sample.config[,2]*colscale
    loc.dists <- sqrt( outer( x, x, "-" )^2 + outer( y, y, "-" )^2 )
    loc.index <- rep(1:nrow(sample.config),sample.config[,4])
    return(loc.dists[loc.index,loc.index])
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
    A <- with( subset(adj, usethese ), Matrix::sparseMatrix( i=i+1L, j=j+1L, x=1.0, dims=c(nrow*ncol,nrow*ncol), symmetric=symmetric ) )
    return(A)
}

###
# stuff to read in and do things with ms output

#' @export nsamples
methods::setGeneric("nsamples", function(x,...) { standardGeneric("nsamples") })

#' @export msseq
methods::setGeneric("msseq", function(x,nloci,...) { standardGeneric("msseq") })

#' @export pimat
methods::setGeneric("pimat", function(x,nloci,...) { standardGeneric("pimat") })

# this is the output of one ms run
methods::setClass("msrun", representation(
        nsamples="numeric",       # number of sampled sequences
        call="character",       # call to ms
        runs="list"             # list of msout objects
    ) )
methods::setMethod("length", signature=c(x="msrun"), definition=function (x) { length(x@runs) } )
methods::setMethod("[[", signature=c(x="msrun",i="ANY"), definition=function (x,i) { x@runs[[i]] } )
methods::setMethod("[[<-", signature=c(x="msrun",i="ANY",value="ANY"), definition=function (x,i,value) { x@runs[[i]]<-value; return(x) } )
methods::setMethod("nsamples", signature=c(x="msrun"), definition=function (x) { x@nsamples } )
methods::setMethod("msseq", signature=c(x="msrun",nloci="numeric"), definition=function (x,nloci) { x@runs <- lapply(x@runs,msseq,nloci=nloci); return(x) } )
methods::setMethod("pimat", signature=c(x="msrun"), definition=function (x) { 
        pimat.seq( do.call( Biostrings::xscat, lapply( x@runs, slot, "sequence" ) ) )
    } )

# this is the output of one ms replicate
methods::setClass("msrep", representation(
        segsites="numeric",     # number of segregating sites
        positions="numeric",    # locations (floats) of the sites
        alleles="matrix"        # binary matrix of alleles
    ) )
methods::setMethod("nsamples", signature=c(x="msrep"), definition=function (x) { nrow(x@alleles) } )
methods::setMethod("msseq", signature=c(x="msrep",nloci="numeric"), definition=function (x,nloci) { make_msseq(x,nloci) } )
methods::setMethod("names", signature=c(x="msrep"), definition=function (x) { rownames(x@alleles) } )
methods::setMethod("pimat", signature=c(x="msrep"), definition=function (x) { pimat.seq(x@sequence) } )

# and, this is what we get after converting to sequence
methods::setClass("msseq", contains="msrep",
    representation(
        sequence="DNAStringSet"     # from Bioconductor::Biostrings
    ) )
methods::setMethod("show",signature=c(object="msseq"), definition=function (object) { show(object@sequence) })

#' @export
to_fasta <- function (msrun,file,...) {
    # 'file' should be a path with a "%" in where the index of the replicate is inserted
    fnames <- sapply( seq_along(msrun@runs), function (k) { gsub("%",k,file) } )
    for (k in seq_along(msrun@runs)) {
        msseq_to_fasta(msrun@runs[[k]],file=fnames[k],...)
    }
    return(fnames)
}

msseq_to_fasta <- function (msseq,file,append=FALSE) {
    dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
    Biostrings::writeXStringSet(msseq@sequence,file=file,append=append)
}

split_fasta <- function (fasta,output=file.path(gsub("[.]fa$","",fasta),"%.fa")) {
    for (fk in seq_along(fasta)) {
        dir.create(dirname(output[fk]),recursive=TRUE,showWarnings=FALSE)
        fseqs <- Biostrings::readDNAStringSet(fasta[fk])
        if (is.null(names(fseqs))) { names(fseqs) <- seq_along(fseqs) }
        for (k in seq_along(fseqs)) {
            outfile <- gsub("%",names(fseqs)[k],output[fk])
            Biostrings::writeXStringSet(fseqs[k],file=outfile)
        }
    }
    return(dirname(output))
}

#' @export
sample_reads <- function (
        infiles,
        coverage,
        error.prob,
        outdirs=gsub("[.]fa$","",infiles),
        count.files=outdirs
    ) {
    # first split up the fasta files into separate ones per individual (sequencer silently uses only the first line)
    # then, in that directory:
    #  use 'sequencer' to sample short reads from input fasta files, writing them out to the files in outfiles;
    #  then map them back to the consensus sequence with 'bwa'
    #  and sort, index, etc with 'samtools'
    #  and count up alleles with 'angsd'
    sequencer.call <- paste("sequencer", "-c", coverage, "-E", error.prob )
    bwa.index.call.head <- paste("bwa index" )
    bwa.call.head <- paste("bwa mem" )
    bwa.call.tail <- paste("|", "samtools view -b - | samtools sort -O bam -T $(mktemp -u tmp.XXXXXX) >")
    samtools.call <- paste("samtools index")
    angsd.call.head <- paste('bash -c "angsd -bam <(echo ')
    angsd.call.mid <- paste(') -out')
    angsd.call.tail <- paste('-minQ 0 -doCounts 1 -dumpCounts 4"')
    # split into indviduals
    for (j in seq_along(infiles)) {
        split_fasta(infiles[j],file.path(outdirs[j],"%.fa"))
        fasta.files <- list.files(outdirs[j],pattern="*[.]fa$",full.names=TRUE)
        reference.file <- fasta.files[1]
        read.files <- gsub("[.]fa$",".reads.fa",fasta.files)
        bam.files <- gsub("[.]fa$",".bam",fasta.files)
        browser()
        for (k in seq_along(fasta.files)) {
            scall <- paste(sequencer.call, infiles[k], ">", read.files[k])
            bicall <- paste(bwa.index.call.head, reference.file)
            bcall <- paste(bwa.call.head, reference.file, read.files[k], bwa.call.tail, bam.files[k])
            samcall <- paste(samtools.call, bam.files[k])
            for (call in c(scall,bicall,bcall,samcall)) {
                cat(call,"\n")
                system(call,ignore.stderr=TRUE)
            }
        }
        acall <- paste( angsd.call.head, paste(bam.files,collapse=" "), angsd.call.mid, count.files[j], angsd.call.tail )
    }
    return(count.files)
}

#' @export
read_msrun <- function (ms.output,
        text=paste(scan(file.path(ms.output,"msoutput.txt"),what='char',sep='\n'),collapse='\n')
    ) {
    msout <- strsplit(text,"\n//\n",fixed=TRUE)[[1]]
    call <- msout[[1]]
    runs <- lapply(msout[-1],read_msrep)
    nsamples <- sapply(runs,nsamples)
    if (length(unique(nsamples))>1) { warning("Not all replicates have the same number of samples?!?") }
    return( new("msrun",
            nsamples=nsamples[1],
            call=call,
            runs=runs
        ) )
}

read_msrep <- function (text) {
    text <- strsplit(text,"\n")[[1]]
    if (!substring(text[1],1,10)=="segsites: ") {
        stop("Expecting 'segsites: ', got", text[1])
    } else {
        segsites <- as.integer(substring(text[1],11))
    }
    if (!substring(text[2],1,11)=="positions: ") {
        stop("Expecting 'positions: ', got", text[2])
    } else {
        positions <- sapply(strsplit(substring(text[2],12)," ")[[1]],as.numeric)
    }
    alleles <- do.call( rbind, strsplit( text[-(1:2)], "" ) )
    if (!all(alleles %in% c("0","1"))) { stop("Found alleles that aren't 0 or 1??") }
    alleles <- ( alleles == "1" )
    # give samples names!
    rownames(alleles) <- paste("s",1:nrow(alleles),sep='')
    return( new( "msrep", segsites=segsites, positions=positions, alleles=alleles ) )
}

make_msseq <- function (msrep,nloci) {
    # Turn a msout object into a msseq object, by picking alleles and adding errors.
    # Moves two mutations occurring in the same bin to the next unmutated bin.
    bases <- c("A","C","G","T")
    refseq <- DNAString( paste(sample(bases,nloci,replace=TRUE),collapse='') )
    loci <- 1+floor( msrep@positions*(nloci-1) )
    if (length(loci)>nloci) { stop("Too many loci!") }
    nn <- 0
    while( nn<100 && any(diff(loci)==0) ) {
        loci[ c(FALSE,diff(loci)==0) ] <- 1+loci[ c(FALSE,diff(loci)==0) ]
    }
    if (nn==100) { stop("Too many mutated sites!") }
    use.loci <- (loci <= nloci)
    if (!all(use.loci)) {
        warning("Some loci fell off the end of the sequence, omitting", sum(!use.loci), "of them.") 
        loci <- loci[use.loci]
    }
    altseq <- DNAString( paste( bases[ 1+((match(strsplit(as.character(refseq[loci]),"")[[1]],bases)+sample.int(3,length(loci),replace=TRUE))%%4) ], collapse="") )
    sequence <- DNAStringSet( list(refseq)[rep(1,nsamples(msrep))] )
    names(sequence) <- names(msrep)
    for (ind in seq_len(nsamples(msrep))) {
        # substitute in alternative alleles
        repthese <- msrep@alleles[ind,use.loci]
        sequence[[ind]][loci[repthese]] <- altseq[repthese]
    }
    return( new( "msseq",
            segsites=msrep@segsites,
            positions=msrep@positions,
            alleles=msrep@alleles,
            sequence=sequence
        ) )
}

add_error_seq <- function (seq,error.prob) {
    # add independent errors to a sequence
    errors <- ( rbinom( nloci, size=1, prob=error.prob ) > 0 )
    nerrors <- sum(errors)
    errseq <- DNAString( paste( bases[ 1+((match(strsplit(as.character(seq[errors]),"")[[1]],bases)+sample.int(3,nerrors,replace=TRUE))%%4) ], collapse="") )
    seq[errors] <- errseq
    return(seq)
}

#' @export
pimat.seq <- function (sequence) {
    # compute the matrix of mean pairwise distances
    Biostrings::stringDist(sequence,method="hamming")/nchar(sequence[[1]])
}
