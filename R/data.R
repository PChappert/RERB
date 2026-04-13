#' AI2B lab COVID-19 Abs dataset
#'
#' A dataset containing all published IgVH and IgVL sequences and linked functional data from the AI2B lab along with associated metadata which includes a codebook describing variables in data.
#'
#' @format A named list with two elements:
#' \describe{
#'   \item{data}{A data frame with 13822 rows and 56 variables. See `metadata$codebook` for full variable descriptions.}
#'   \item{metadata}{A named list containing:
#'     \describe{
#'       \item{donors}{Data frame (32 × 22) with donor-level information}
#'       \item{sampling}{Data frame (32 × 16) with sampling time points}
#'       \item{codebook}{Data frame (56 × 3) describing variables in `data`, `metadata$donors` and `metadata$sampling`}
#'     }
#'   }
#' }
#' 
#' @source currently includes sequences from:
#' Sokal_Cell_2021 (PMID: 33571429)
#' Sokal_Immunity_2021 (PMID: 35483354)
#' Sokal_Immunity_2022 (PMID: 35483354)
#' Sokal_JEM_2022 (PMID: 36342455)
#' Sokal_Immunity_2023 (PMID: 37543032)
#' @examples
#' data(AI2B_COVID_dataset)
#' AI2B_COVID_dataset$data
#' AI2B_COVID_dataset$metadata$donors
#' AI2B_COVID_dataset$metadata$sampling
#' AI2B_COVID_dataset$metadata$codebook
"AI2B_COVID_dataset"