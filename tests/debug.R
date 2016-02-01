library(Matrix)
library(msarg)

# Was having issues with 
#   " infinite time to next event. Negative growth rate in last time interval or non-communicating subpops."
# in this situation, because of zeros on the diagonal of @M.

M1 <- new("dgCMatrix"
            , i = c(1L, 3L, 4L, 0L, 2L, 3L, 4L, 5L, 1L, 4L, 5L, 0L, 1L, 4L, 6L, 7L, 
                    0L, 1L, 2L, 3L, 5L, 6L, 7L, 8L, 1L, 2L, 4L, 7L, 8L, 3L, 4L, 7L, 
                    3L, 4L, 5L, 6L, 8L, 4L, 5L, 7L)
            , p = c(0L, 3L, 8L, 11L, 16L, 24L, 29L, 32L, 37L, 40L)
            , Dim = c(9L, 9L)
            , Dimnames = list(NULL, NULL)
            , x = rep(7e-5, 40)
            , factors = list()
        )

M2 <- new("dgCMatrix"
        , i = c(0L, 1L, 3L, 4L, 0L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 4L, 5L, 0L, 
                1L, 3L, 4L, 6L, 7L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 
                4L, 5L, 7L, 8L, 3L, 4L, 6L, 7L, 3L, 4L, 5L, 6L, 7L, 8L, 4L, 5L, 
                7L, 8L)
        , p = c(0L, 4L, 10L, 14L, 20L, 29L, 35L, 39L, 45L, 49L)
        , Dim = c(9L, 9L)
        , Dimnames = list(NULL, NULL)
        , x = c(0, 7e-05, 7e-05, 6e-06, 7e-05, 
                0, 7e-05, 6e-06, 7e-05, 6e-06, 7e-05, 
                0, 6e-06, 7e-05, 7e-05, 6e-06, 
                0, 7e-05, 7e-05, 6e-06, 6e-06, 7e-05, 6e-06, 7e-05, 
                0, 7e-05, 6e-06, 7e-05, 6e-06, 6e-06, 7e-05, 7e-05, 
                0, 6e-06, 7e-05, 7e-05, 6e-06, 
                0, 7e-05, 6e-06, 7e-05, 6e-06, 7e-05, 
                0, 7e-05, 6e-06, 7e-05, 7e-05, 0)
        , factors = list()
    )

dem <- new( "demography",
        popStates = list(
                new("gridArray"
                    , npop = c(3L, 3L, 1L)
                    , N = rep(3000,9)
                    , G = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
                    , M = M1
                ),
                new("gridArray"
                    , npop = c(3L, 3L, 1L)
                    , N = rep(3000,9)
                    , G = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
                    , M = M2
                ) ),
        t=20000 )

sample.config <-  structure(c(3, 3, 1, 2, 1, 1, 7, 3), .Dim = c(2L, 4L), .Dimnames = list(
                    NULL, c("row", "col", "layer", "n")))

run_ms( dem, nsamp=sample.config, trees=TRUE, outdir="ms_debug", nreps=3 )
