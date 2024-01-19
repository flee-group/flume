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

#' Algal metacommunity N:P dataset
#'
#' @format A list
#' \describe{
#'	\item{niche_params}{A data.frame giving the niche location, breadth, and scale (in terms of N:P) for 40 species}
#'	\item{r0}{Starting resource state for N and P}
#'	\item{sp0}{Starting locations for all species}
#' }
"algae"

#' Flume river network and associated data for the Kamp River catchment in Austria
#'
#' @format A `list`, with members:
#' \describe{
#'   \item{network}{A [river_network]}
#' }
"kamp"

#' Flume river networks for 5 different rivers
#'
#' @format A `list`, with members:
#' \describe{
#'   \item{kamp}{Kamp river, Austria}
#'   \item{mara}{Mara river, Kenya}
#'   \item{thur}{Thur river, Switzerland}
#'   \item{vjosa}{Vjosa river, Greece & Albania}
#'   \item{ybbs}{Ybbs river, Austria}
#' }
"flume_networks"
