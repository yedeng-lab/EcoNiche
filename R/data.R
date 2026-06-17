#' Example plant community abundance table
#'
#' A forest plant community dataset for demonstrating the EcoNiche workflow.
#' Rows are plant species and columns are samples. Values are non-negative
#' abundance counts.
#'
#' @format A numeric matrix with 1038 plant species in rows and 290 samples in
#' columns.
#' @seealso \code{\link{plant_env}}, \code{\link{plant_group}}
"plant_otu"

#' Example plant environmental table
#'
#' Environmental variables for the samples in \code{\link{plant_otu}}. Rows are
#' samples and columns are environmental variables.
#'
#' @format A data frame with 290 rows and 11 variables:
#' \describe{
#'   \item{TN}{Total nitrogen.}
#'   \item{TK}{Total potassium.}
#'   \item{TP}{Total phosphorus.}
#'   \item{SOC}{Soil organic carbon.}
#'   \item{AP}{Available phosphorus.}
#'   \item{AN}{Available nitrogen.}
#'   \item{Moisture}{Soil moisture.}
#'   \item{AMT}{Annual mean temperature.}
#'   \item{AMP}{Annual mean precipitation.}
#'   \item{pH}{Soil pH.}
#'   \item{abspH}{Absolute deviation from neutral pH.}
#' }
#' @seealso \code{\link{plant_otu}}, \code{\link{plant_group}}
"plant_env"

#' Example plant sample grouping
#'
#' Group labels for the samples in \code{\link{plant_otu}} and
#' \code{\link{plant_env}}. The vector names are sample IDs.
#'
#' @format A named character vector with one group label per sample.
#' @seealso \code{\link{plant_otu}}, \code{\link{plant_env}}
"plant_group"
