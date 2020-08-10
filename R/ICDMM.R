#' ICDMM: Deterministic malaria transmission model
#'
#' ICDMM implements the Imperial College malaria transmission model
#' using the R package odin. A related suite of tools for setting up
#' different regional settings and model types is provided provides a flexible grammar of data manipulation.
#'
#' @useDynLib ICDMM
#' @importFrom stats rlnorm
#' @importFrom utils adist
#' @importFrom ggplot2 ggplot
#' @importFrom rlang .data
#' @importFrom odin odin
#' @importFrom dde difeq
NULL

cache <- new.env()
