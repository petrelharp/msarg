
context("refugia")

set.seed(42)

grid.nrow <- 2                 # width of grid
grid.ncol <- 20                # height of grid
glacier.width <- 17       # with of the glacier (in the middle)
mig.rate <- 1             # migration rate between neighbors in the grid
growth.time <- 200          # units of time (in Ne) it took to re-expand post-glacier
glacier.end <- 500          # units of time ago (in Ne) that glaciation ended
glacier.begin <- 1000        # units of time ago (in Ne) that glaciation began
pop.density <- 1+rpois( grid.nrow*grid.ncol, 10 )

dem <- ms_demog( grid_array(nlayers=1,nrow=grid.nrow,ncol=grid.ncol,
                N=pop.density, mig.rate=mig.rate) )

mask <- matrix(1,nrow=nrow(dem),ncol=ncol(dem))
mask[ abs(col(mask)-ncol(dem)/2) < glacier.width/2-0.1] <- 0
dem <- add_to_demography( dem, tnew=glacier.end, fn=modify_grid_layer, layer=1, dN=mask )

speed <- 1.2 * ( glacier.width ) / growth.time / 10
width <- 4
dem <- logistic_interpolation( dem, 
    t.end=glacier.end-growth.time,   # time ago expansion ended
    t.begin=glacier.end,             # time ago exapnsion began
    nsteps=10,                        # number of time steps to pass to ms
    speed=speed, width=width )       # parameters; could also pass sigma= and r= .

dem <- add_to_demography( dem, tnew=glacier.begin, pop=1 )

sample.config <- sample_locations( dem, n=30, each=2 )

ms.output <- "ms_testing"
test_that( "running with reexpansion", {
    expect_error( run_ms( dem, nsamp=sample.config, trees=TRUE, nreps=3, outdir=ms.output ), NA );
    expect_equal_to_reference( trees_from_ms( file.path(ms.output,"msoutput.txt") ), "saved-refugia-trees.Rd" )
    } )

unlink(file.path(ms.output,"msoutput.txt"))
unlink(file.path(ms.output,"msarg.txt"))
unlink(ms.output)

if (FALSE) {  # takes too long, not really testing anything

    context("larger pop sizes")

    big.dem <- dem
    for (k in 1:length(big.dem)) { big.dem[[k]]@N <- big.dem[[k]]@N * 1e4 }

    test_that( "reexpansion, big pops", {
        expect_that( run_ms( big.dem, nsamp=sample.config, trees=TRUE, nreps=3, outdir=ms.output ),
                     not(throws_error()) );
        expect_equal_to_reference( trees_from_ms( file.path(ms.output,"msoutput.txt") ), "saved-refugia-trees-2.Rd" )
        } )


    unlink(file.path(ms.output,"msoutput.txt"))
    unlink(file.path(ms.output,"msarg.txt"))
    unlink(ms.output)



    context("smaller migration rates")

    slow.dem <- dem
    for (k in 1:length(slow.dem)) { slow.dem[[k]]@M@x <- slow.dem[[k]]@M@x * 1e-8 }

    test_that( "reexpansion, small migration rates", {
        expect_that( 
                run_ms( slow.dem, nsamp=sample.config, trees=TRUE, nreps=3, outdir=ms.output )
                , not(throws_error()) );
        expect_equal_to_reference( 
                trees_from_ms( file.path(ms.output,"msoutput.txt") )
                , "saved-refugia-trees-3.Rd" )
        } )


    unlink(file.path(ms.output,"msoutput.txt"))
    unlink(file.path(ms.output,"msarg.txt"))
    unlink(ms.output)

}
