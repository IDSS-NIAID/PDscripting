# utils.R
# Utility functions for the package

#' read_sample_node_args
#'
#' Read in a sample node arguments file provided with the package
#'
#' @param set Numeric value, indicating the sample node_args file to read in
#'
#' @return A list containing the sample node arguments
#' @export
#' @importFrom rjson fromJSON
#' @importFrom stringr str_replace
read_sample_node_args <- function(set = NULL)
{
  # check that the set is valid
  file <- system.file('extdata', paste0('node_args', set, '.json'), package = "PDscripting")

  if(file == '')
  {
    sets <- system.file("extdata", package = "PDscripting") |>
      list.files() |>
      grep(pattern = "node_args", value = TRUE) |>
      str_replace(pattern = "node_args", replacement = "") |>
      str_replace(pattern = ".json", replacement = "")

    stop("Please provide a valid set number from the following: ", paste(sets, collapse = ", "))
  }

  # read in the file
  node_args <- rjson::fromJSON(file = file)
  node_args$Tables$DataFile <- eval(parse(text = node_args$Tables$DataFile))

  return(node_args)
}
