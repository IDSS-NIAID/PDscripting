# PR_fitler
# Filter replicates for valid protein ratios

# documentation on custom script integration: https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1325195659v1


#### Setup ####
#' Requested tables and columns: "Proteins: Accession, Abundances Normalized, Found in Samples"

debug_dir <- 'C:/Users/cormierkw/Downloads'


library(rjson)
library(dplyr)
library(readr)


#### Read in node args from PD ####

node_args <- fromJSON(file = commandArgs(trailingOnly = TRUE)[1])

if(FALSE)
{
  # when running manually, navigate to the R directory if we aren't already there
  root <- here::here()
  file.path(root, 'R') |>
    setwd()
  
  node_args <- fromJSON(file = 'data/UC1/node_args.json')
}


#### Fetch and wrangle data ####

dat <- read_delim(node_args$Tables[[1]]$DataFile)

write_delim(dat, file = file.path(debug_dir, 'to_pd.tsv'))


#### Write output to PD ####

list(CurrentWorkflowID = node_args$CurrentWorkflowID,
     Tables = list(list(TableName = "Proteins",
                        DataFile  = file.path(debug_dir, 'to_pd.tsv'),
                        DataFormat = "CSV",
                        Options = '{}',
                        ColumnDescriptions = list(list(ColumnName = "Accession",
                                                       ID = "ID",
                                                       DataType = "String",
                                                       Options = "{}"),
                                                  list(ColumnName = "AbundancesFiltered",
                                                       ID = "",
                                                       DataType = "Float",
                                                       Options = "")
                                                 )
                       )
                  )
     ) |>
  toJSON(indent = 2) |>
  cat(file = node_args$ExpectedResponsePath)
