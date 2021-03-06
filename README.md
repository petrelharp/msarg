**msarg** : some `R` infrastructure to run `ms` simulations on grids
====================================================================

Here are some S4 classes and some methods for building up demographic models on spatial grids
that change over time,
along with methods to write these out as command-line parameters for [`ms`](http://home.uchicago.edu/rhudson1/source/mksamples.html).

Have a look at the examples:

1. Appearance, and disappearance, of [a barrier](http://petrelharp.github.io/msarg/barrier-howto.html)
2. [Parallel universes](parallel-universes-howto.html)
3. Sudden appearance, and gradual retreat, of [a glacier](http://petrelharp.github.io/msarg/glaciation-howto.html)

... and there's a bit more documentation, sort of in [this document](http://petrelharp.github.io/msarg/using-msarg.html).

intalling
---------

This should do it:

```
devtools::install_github("petrelharp/msarg")
```

getting ms
----------

Go to [Hudson's site](http://home.uchicago.edu/rhudson1/source/mksamples.html) and download it.  The proper citation for `ms` is: *"Hudson, R. R. (2002) Generating samples under a Wright-Fisher neutral model. Bioinformatics 18:337-8."*
