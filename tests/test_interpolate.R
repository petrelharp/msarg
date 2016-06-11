source("../msarg.R")

dem <- ms_demog( grid_array(nlayers=1,nrow=30,ncol=1,N=1,mig.rate=1) )
mask <- matrix(1,nrow=nrow(dem),ncol=ncol(dem))
mask[ row(mask)>2 ] <- 0
dem <- add_to_demography( dem, tnew=10, fn=modify_grid_layer, layer=1, dN=mask )

speed <- 20 / 10
width <- 3
dem <- logistic_interpolation( dem, 
    t.end=1,   # time ago expansion ended
    t.begin=10,             # time ago exapnsion began
    nsteps=10,                        # number of time steps to pass to ms
    speed=speed, width=width )       # parameters; could also pass sigma= and r= .


plot(dem)

