# wrangling.R
# functions for wrangling data

#' PD_colTypes
#' 
#' Fetch column types from PD node_args json file
#' 
#' @param colDescriptions list of column descriptions from PD node_args
#' 
#' @description
#' This function reads in all column names from the node_args.json file and returns a tibble with the following columns:
#' - name: original column name
#' - type: expected column type
#' - new_name: new column name
#' 
#' @return tibble with column names and types
#' @export
#' @importFrom dplyr case_when
#' @importFrom purrr map_chr
#' @importFrom stringr str_replace str_replace_all
PD_colTypes <- function(colDescriptions)
{
  tibble(name = map_chr(colDescriptions, ~ .x$ColumnName),
         type = map_chr(colDescriptions, ~ 
                          case_when(.x$ID       == "ID"      ~ 'c', # IDs are large negative integers
                                    .x$DataType == "Float"   ~ 'n',
                                    .x$DataType == "String"  ~ 'c',
                                    .x$DataType == "Int"     ~ 'i',
                                    TRUE                     ~ '?')),
         new_name = case_when(grepl("Abundances Normalized", name) ~ str_replace(name, ' ', '.') |> # convert to "Abundances.Normalized <sample Id>.Sample"
                                str_replace(' ', '_') |>       # convert to "Abundances.Normalized_<sample Id>.Sample"
                                str_replace_all(' ', '.'),     # convert to "Abundances.Normalized_<sample.Id>.Sample"
                              
                              grepl("Found in Sample",       name) ~ str_replace(name, ' ', '.') |> # convert to "Found.in.Samples <sample Id>.Sample"
                                str_replace(' ', '.') |> 
                                str_replace(' ', '_') |>       # convert to "Found.in.Samples_<sample Id>.Sample"
                                str_replace_all(' ', '.'),     # convert to "Found.in.Samples_<sample.Id>.Sample"
                              
                              grepl("Abundance Ratio",       name) ~ str_replace(name, ' ', '.') |> # convert to "Abundance.Ratio* <sample Id>  <sample Id>"
                                str_replace('Ratio log2', 'Ratio.log2') |>
                                str_replace('Ratio Adj', 'Ratio.Adj') |>
                                str_replace(' P-Value', '.P-Value') |>
                                str_replace(' ', '_') |>       # convert to "Abundance.Ratio*_<sample Id>  <sample Id>"
                                str_replace('  ', '/') |>      # convert to "Abundance.Ratio*_<sample Id>/<sample Id>"
                                str_replace_all(' ', '.'),     # convert to "Abundance.Ratio*_<sample.Id>/<sample.Id>"
                              
                              grepl('Abundances Grouped',    name) ~ str_replace(name, ' ', '.') |> # convert to "Abundances.Grouped* <sample Id>"
                                str_replace('Grouped Count', 'Grouped.Count') |>
                                str_replace(' ', '_') |>       # convert to "Abundances.Grouped*_<sample Id>"
                                str_replace_all(' ', '.'),     # convert to "Abundances.Grouped*_<sample.Id>"
                              
                              TRUE ~ name))
}


#' PD_clean
#' 
#' Clean up data
#' 
#' @param file file to read in from PD and clean
#' @param colTypes column types data.frame from `PD_colTypes()`
#' 
#' @return data frame with cleaned data
#' @export
#' @importFrom dplyr mutate
#' @importFrom readr read_delim
#' @importFrom stringr str_trim
PD_clean <- function(file, colTypes)
{
  # read in data
  dat <- read_delim(file, col_types = colTypes$type,
                    na = c('', 'NA', 'n/a', ' n/a'))
  
  # some columns have extra spaces in the data - double check and remove them
  for(i in 1:ncol(dat))
  {
    # expected character columns
    if(colTypes$type[i] == 'c')
    {
      dat[[i]] <- str_trim(dat[[i]])
    }
    
    # expected numeric columns but read in as character - probably has spaces
    if(colTypes$type[i] == 'n' & is(dat[[i]], 'character'))
    {
      dat[[i]] <- str_trim(dat[[i]]) |> as.numeric()
    }
    
    # expected integer columns but read in as character
    # many of these seem to be more like factors than integers, but not all :P
    if(colTypes$type[i] == 'i' & is(dat[[i]], 'character'))
    {
      dat[[i]] <- str_trim(dat[[i]])                   # remove white space
      dat[[i]] <- ifelse(dat[[i]] == '', NA, dat[[i]]) # remove blanks
      
      # if all values are actually integers, convert to integer
      if(all(!is.na(suppressWarnings(as.integer(names(table(dat[[i]])))))))
      {
        dat[[i]] <- as.integer(dat[[i]])
      }
    }else if(colTypes$type[i] == 'i' & is(dat[[i]], 'numeric')){
      dat[[i]] <- as.integer(dat[[i]])
    }
  }
  
  
  # make column names a little easier to work with
  colnames(dat) <- colTypes$new_name
  
  return(dat)
}


#' write_PD_response
#' 
#' Write PD response file
#' 
#' @param node_args node args from PD
#' @param dat filtered data frame to write out
#' @param colTypes column types data.frame from `PD_colTypes()`
#' 
#' @description
#' This function writes out the node response file to the location defined in `node_args$ExpectedResponsePath`
#' 
#' @return Invisibly returns node_response JSON list
#' @export
#' @importFrom dplyr mutate
#' @importFrom jsonlite toJSON
#' @importFrom stringr str_replace str_replace_all
write_PD_response(node_args, dat, colTypes, startsWith)
{
  # create node response
  node_response = list(CurrentWorkflowID = node_args$CurrentWorkflowID,
                       Tables = list(list(TableName = "Proteins",
                                          DataFile  = file.path(dirname(node_args$Tables[[1]]$DataFile),
                                                                'to_pd.tsv'),
                                          DataFormat = "CSV",
                                          Options = NULL,
                                          ColumnDescriptions = list())
                       )
  )
  
  
  # update colTypes with new names
  for(cols in startsWith)
  {
    colTypes <- mutate(colTypes,
                       new_name = str_replace(new_name, fixed(cols), paste0('Filtered.', cols)))
  }
  
  # check that all column names were updated correctly
  if(!all(colnames(dat_filtered) %in% colTypes$new_name))
  {
    stop('Not all column names in `dat_filtered` were updated correctly')
  }
  
  # loop over new column names
  nextCol <- 1
  for(nn in names(dat_filtered))
  {
    # figure out which old name we are looking for
    on <- colTypes$name[colTypes$new_name == nn]
    
    for(j in 1:length(node_args$Tables[[1]]$ColumnDescriptions))
    {
      # when we find it, copy metadata over and update
      if(node_args$Tables[[1]]$ColumnDescriptions[[j]]$ColumnName == on)
      {
        node_response$Tables[[1]]$ColumnDescriptions[[nextCol]] <- node_args$Tables[[1]]$ColumnDescriptions[[j]]
        
        node_response$Tables[[1]]$ColumnDescriptions[[nextCol]]$ColumnName <- nn
        
        nextCol <- nextCol + 1
      }
    }
  }
  
  # write out node response
  toJSON(node_response, indent = 2) |>
    str_replace_all(fixed('"{}"'), '{}') |> # for some reason it wants to quote empty lists
    str_replace_all(fixed("[\n\n          ]\n"), "{}") |> # for some reason we end up with these, too
    cat(file = node_args$ExpectedResponsePath)
  
  invisible(node_response)
}
