# flume
An R-package for a FLUvial Meta Ecosystem Model

## Installation

First steps:
Flume requires R 4 or later, so please upgrade if you are still on an early version of R.

There is a companion package in development, [watershed](https://github.com/flee-group/watershed), that is useful for importing river networks into flume (also see the flume vignette about [importing a river network](https://flee-group.github.io/flume/network_import.html)).

```r
## requires a recent version of the remotes package
remotes::install_github("flee-group/flume")

## if you are out of date, try:
remotes::install_github("flee-group/flume", ref="main")
```

## Getting started

There are a number of vignettes available to help you get a flume model going. Start with [a simple sumulation](https://flee-group.github.io/flume/simple_sim.html) to see how the model is structured. Then you can get more information on [importing a river network](https://flee-group.github.io/flume/network_import.html) and [designing metacommunities](https://flee-group.github.io/flume/metacommunities.html).
