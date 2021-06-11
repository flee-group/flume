#' flume: A package for modelling FLUvial Meta-Ecosystems
#'
#' Implementation of a theoretical/mechanistic model for exploring fluvial meta-ecosystems.
#'
#' @section Key functions:
#'
#' * [metacommunity()] Create species pools and associated attributes
#' * [river_network()] Build a river network object
#' * [flume()] Creates a model
#' * [run_simulation()] Runs the simulation
#'
#' @docType package
#' @name flume_package
NULL

#' Algal metacommunity N:P example dataset
#'
#' @format A [flume()]
#' \describe{
#'   \item{species}{Five algal ASVs}
#'   \item{sites}{Sampled locations in the Bodingbach catchment}
#' }
"algae"

#' Flume river network and associated data for the Kamp River catchment in Austria
#'
#' @format A `list`, with members:
#' \describe{
#'   \item{network}{A [river_network]}
#' }
"kamp"

