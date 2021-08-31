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
#' @format A list
#' \describe{
#'	\item{niches}{A table giving the niche location and breadth (in terms of N:P) for five species}
#'	\item{adjacency}{Adjacency matrix for the river network}
#'	\item{discharge}{Discharge at each node in the river network}
#'	\item{cs_area}{Cross-sectional area for each node in the network}
#'	\item{layout}{Plotting layout for the river network}
#'	\item{network}{River network object constructed from the algae data}
#'	\item{r0}{Starting resource state for N and P}
#'	\item{sp0}{Starting locations for all species}
#'	\item{metacommunity}{A metacommunity object constructed from all the data}
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
