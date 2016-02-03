library(Matrix)
library(msarg)

context("dealing with zeros correctly")

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

# this one has zeros on the diagonal
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


# this version does not have zeros on the diagonal
M3 <- as( as.matrix( M2 ), "dgCMatrix" )
dem3 <- dem
dem3[[2]]@M <- M3

sample.config <-  structure(c(3, 3, 1, 2, 1, 1, 7, 3), .Dim = c(2L, 4L), .Dimnames = list(
                    NULL, c("row", "col", "layer", "n")))

args1 <- msarg( dem, nsamp=sample.config )
args3 <- msarg( dem3, nsamp=sample.config )

# we do NOT expect this, because the redundancy-checking isn't that good, adds '-em' flags it doesn't need to (but no harm)
# test_that( "same matrices different structure",
#           expect_identical( args1, args3 )
#       )


run1 <- run_ms( dem, nsamp=sample.config, trees=TRUE, tofile=FALSE, nreps=3, seed=c(1,2,3) )
run3 <- run_ms( dem3, nsamp=sample.config, trees=TRUE, tofile=FALSE, nreps=3, seed=c(1,2,3) )

test_that( "simulation: 1",
        expect_equal(
                                  paste(run1, collapse="; ")
                          ,
"ms 10 3 -T -seed 1 2 3 -I 9 0 0 7 0 0 3 0 0 0 -n 1 3000 -n 2 3000 -n 3 3000 -n 4 3000 -n 5 3000 -n 6 3000 -n 7 3000 -n 8 3000 -n 9 3000 -g 1 0 -g 2 0 -g 3 0 -g 4 0 -g 5 0 -g 6 0 -g 7 0 -g 8 0 -g 9 0 -m 2 1 7e-05 -m 4 1 7e-05 -m 5 1 7e-05 -m 1 2 7e-05 -m 3 2 7e-05 -m 4 2 7e-05 -m 5 2 7e-05 -m 6 2 7e-05 -m 2 3 7e-05 -m 5 3 7e-05 -m 6 3 7e-05 -m 1 4 7e-05 -m 2 4 7e-05 -m 5 4 7e-05 -m 7 4 7e-05 -m 8 4 7e-05 -m 1 5 7e-05 -m 2 5 7e-05 -m 3 5 7e-05 -m 4 5 7e-05 -m 6 5 7e-05 -m 7 5 7e-05 -m 8 5 7e-05 -m 9 5 7e-05 -m 2 6 7e-05 -m 3 6 7e-05 -m 5 6 7e-05 -m 8 6 7e-05 -m 9 6 7e-05 -m 4 7 7e-05 -m 5 7 7e-05 -m 8 7 7e-05 -m 4 8 7e-05 -m 5 8 7e-05 -m 6 8 7e-05 -m 7 8 7e-05 -m 9 8 7e-05 -m 5 9 7e-05 -m 6 9 7e-05 -m 8 9 7e-05 -em 20000 2 1 7e-05 -em 20000 4 1 7e-05 -em 20000 5 1 6e-06 -em 20000 1 2 7e-05 -em 20000 3 2 7e-05 -em 20000 4 2 6e-06 -em 20000 5 2 7e-05 -em 20000 6 2 6e-06 -em 20000 2 3 7e-05 -em 20000 5 3 6e-06 -em 20000 6 3 7e-05 -em 20000 1 4 7e-05 -em 20000 2 4 6e-06 -em 20000 5 4 7e-05 -em 20000 7 4 7e-05 -em 20000 8 4 6e-06 -em 20000 1 5 6e-06 -em 20000 2 5 7e-05 -em 20000 3 5 6e-06 -em 20000 4 5 7e-05 -em 20000 6 5 7e-05 -em 20000 7 5 6e-06 -em 20000 8 5 7e-05 -em 20000 9 5 6e-06 -em 20000 2 6 6e-06 -em 20000 3 6 7e-05 -em 20000 5 6 7e-05 -em 20000 8 6 6e-06 -em 20000 9 6 7e-05 -em 20000 4 7 7e-05 -em 20000 5 7 6e-06 -em 20000 8 7 7e-05 -em 20000 4 8 6e-06 -em 20000 5 8 7e-05 -em 20000 6 8 6e-06 -em 20000 7 8 7e-05 -em 20000 9 8 7e-05 -em 20000 5 9 6e-06 -em 20000 6 9 7e-05 -em 20000 8 9 7e-05 ; 1 2 3; ; //; ((10:934.296,(8:237.254,9:237.254):697.042):2572.808,(6:2034.242,(2:1441.670,((4:96.124,(1:95.365,3:95.365):0.759):162.988,(5:127.675,7:127.675):131.437):1182.558):592.572):1472.861);; ; //; (9:28284.086,(((7:94.415,(1:20.146,2:20.146):74.269):790.948,(3:260.163,6:260.163):625.200):11807.425,((8:115.966,10:115.966):8551.475,(4:225.016,5:225.016):8442.425):4025.347):15591.299);; ; //; ((10:2151.344,(1:1170.637,((2:196.277,(4:146.457,5:146.457):49.820):181.111,(7:326.394,(3:177.582,6:177.582):148.813):50.994):793.249):980.707):229937.109,(8:26026.879,9:26026.879):206061.578);"
        ) )

test_that( "simulation: 2",
        expect_equal(
                                  paste(run3, collapse="; ")
                          ,
"ms 10 3 -T -seed 1 2 3 -I 9 0 0 7 0 0 3 0 0 0 -n 1 3000 -n 2 3000 -n 3 3000 -n 4 3000 -n 5 3000 -n 6 3000 -n 7 3000 -n 8 3000 -n 9 3000 -g 1 0 -g 2 0 -g 3 0 -g 4 0 -g 5 0 -g 6 0 -g 7 0 -g 8 0 -g 9 0 -m 2 1 7e-05 -m 4 1 7e-05 -m 5 1 7e-05 -m 1 2 7e-05 -m 3 2 7e-05 -m 4 2 7e-05 -m 5 2 7e-05 -m 6 2 7e-05 -m 2 3 7e-05 -m 5 3 7e-05 -m 6 3 7e-05 -m 1 4 7e-05 -m 2 4 7e-05 -m 5 4 7e-05 -m 7 4 7e-05 -m 8 4 7e-05 -m 1 5 7e-05 -m 2 5 7e-05 -m 3 5 7e-05 -m 4 5 7e-05 -m 6 5 7e-05 -m 7 5 7e-05 -m 8 5 7e-05 -m 9 5 7e-05 -m 2 6 7e-05 -m 3 6 7e-05 -m 5 6 7e-05 -m 8 6 7e-05 -m 9 6 7e-05 -m 4 7 7e-05 -m 5 7 7e-05 -m 8 7 7e-05 -m 4 8 7e-05 -m 5 8 7e-05 -m 6 8 7e-05 -m 7 8 7e-05 -m 9 8 7e-05 -m 5 9 7e-05 -m 6 9 7e-05 -m 8 9 7e-05 -em 20000 5 1 6e-06 -em 20000 4 2 6e-06 -em 20000 6 2 6e-06 -em 20000 5 3 6e-06 -em 20000 2 4 6e-06 -em 20000 8 4 6e-06 -em 20000 1 5 6e-06 -em 20000 3 5 6e-06 -em 20000 7 5 6e-06 -em 20000 9 5 6e-06 -em 20000 2 6 6e-06 -em 20000 8 6 6e-06 -em 20000 5 7 6e-06 -em 20000 4 8 6e-06 -em 20000 6 8 6e-06 -em 20000 5 9 6e-06 ; 1 2 3; ; //; ((10:934.296,(8:237.254,9:237.254):697.042):2572.808,(6:2034.242,(2:1441.670,((4:96.124,(1:95.365,3:95.365):0.759):162.988,(5:127.675,7:127.675):131.437):1182.558):592.572):1472.861);; ; //; (9:69383.031,(((7:94.415,(1:20.146,2:20.146):74.269):790.948,(3:260.163,6:260.163):625.200):11807.425,((8:115.966,10:115.966):8551.475,(4:225.016,5:225.016):8442.425):4025.347):56690.242);; ; //; ((4:3045.721,(1:705.171,(3:210.099,7:210.099):495.072):2340.550):32883.367,((6:112.612,(2:18.653,5:18.653):93.959):28385.537,(8:1534.348,(9:4.857,10:4.857):1529.491):26963.801):7430.938);"
        ) )

