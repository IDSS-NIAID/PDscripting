# PR_fitler
# Filter replicates for valid protein ratios

# documentation on custom script integration: https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1325195659v1


#### Setup ####
#' Requested tables and columns: "Proteins: Accession, Abundances Normalized, Found in Samples"


library(rjson)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)

# proportion of samples with non-missing, above-detection-limit values required to keep abundance measures
prop_good <- 0.5

# detection limit
llod <- 1e6


#### Read in node args from PD ####

node_args <- fromJSON(file = commandArgs(trailingOnly = TRUE)[1])

if(FALSE)
{
  # for manual debugging
  node_args <- fromJSON(file = 'debug/node_args.json')
  
  # change paths to `debug/`
  node_args$Tables[[1]]$DataFile <- str_replace(node_args$Tables[[1]]$DataFile, 
                                                fixed(dirname(node_args$Tables[[1]]$DataFile)),
                                                'debug')
  node_args$ExpectedResponsePath <- str_replace(node_args$ExpectedResponsePath, 
                                                fixed(dirname(node_args$ExpectedResponsePath)),
                                                'debug')
}


#### Fetch and wrangle data ####

# fetch expected column types
colTypes <- tibble(name = map_chr(node_args$Tables[[1]]$ColumnDescriptions, ~ .x$ColumnName),
                   type = map_chr(node_args$Tables[[1]]$ColumnDescriptions, ~ 
                                    case_when(.x$ID       == "ID"      ~ 'c', # IDs are large negative integers
                                              .x$DataType == "Float"   ~ 'n',
                                              .x$DataType == "String"  ~ 'c',
                                              .x$DataType == "Int"     ~ 'i',
                                              TRUE                     ~ '?')),
                   new_name = case_when(grepl("Abundances Normalized", name) ~ str_replace(name, ' ', '.') |> # convert to "Abundances.Normalized <sample Id>"
                                                                               str_replace(' ', '_') |>       # convert to "Abundances_Normalized_<sample Id>"
                                                                               str_replace_all(' ', '.'),     # convert to "Abundances.Normalized_<sample.Id>"
                                         grepl("Found in Sample", name)      ~ str_replace(name, ' ', '.') |> # convert to "Found.in.Samples <sample Id>"
                                                                               str_replace(' ', '.') |> 
                                                                               str_replace(' ', '_') |>       # convert to "Found.in.Samples_<sample Id>"
                                                                               str_replace_all(' ', '.'),     # convert to "Found.in.Samples_<sample.Id>"
                                        TRUE ~ name))


# read in data
dat <- read_delim(node_args$Tables[[1]]$DataFile,
                  col_types = colTypes$type,
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
  # these seem to be more like factors than integers - leave as characters
  if(colTypes$type[i] == 'i' & is(dat[[i]], 'character'))
  {
    dat[[i]] <- str_trim(dat[[i]])
  }
}


# make column names a little easier to work with
colnames(dat) <- colTypes$new_name


# convert to long format for summarization
dat_long <- dat |>
  pivot_longer(cols = starts_with("Abundances") | starts_with("Found"),
               names_to = c(".value", "sample_replicate"),
               names_pattern = "(.*)_(.*)") |>
  mutate(sample = str_replace(sample_replicate, 'F\\d+\\.Sample\\.', '') |>
                  str_replace('BR\\d+\\.', ''))


# filter out samples with too many missing or low values
dat_filtered <- dat_long |>

  # count good abundances per sample
  group_by(`Proteins Unique Sequence ID`, Accession, sample) |>
  mutate(thresh_good = ceiling(n() * prop_good),
         ngood = sum(!is.na(Abundances.Normalized) & Abundances.Normalized > llod)) |>
  ungroup() |>
  
  # filter any samples that don't meet the threshold
  mutate(Abundances.Normalized = ifelse(ngood >= thresh_good, Abundances.Normalized, NA)) |>
  
  # pivot back to wide format
  pivot_wider(names_from = sample_replicate,
              values_from = c("Abundances.Normalized", "Found.in.Sample")) |>
  
  # drop unnecessary columns
  select(-Accession, -ngood, -thresh_good, -sample, -starts_with("Found"))


# update column names and write file
colnames(dat_filtered) <- str_replace(colnames(dat_filtered), 'Normalized', 'Normalized.Filtered')

out_dir <- dirname(node_args$Tables[[1]]$DataFile)
write_delim(dat_filtered, file = file.path(out_dir, 'to_pd.tsv'))


#### Write output to PD ####

node_response = list(CurrentWorkflowID = node_args$CurrentWorkflowID,
                     Tables = list(list(TableName = "Proteins",
                                        DataFile  = file.path(out_dir, 'to_pd.tsv'),
                                        DataFormat = "CSV",
                                        Options = '{}',
                                        ColumnDescriptions = list())
                                  )
                    )

nextCol <- 1
for(i in 1:length(node_args$Tables[[1]]$ColumnDescriptions))
{
  if(node_args$Tables[[1]]$ColumnDescriptions[[i]]$ColumnName == "Proteins Unique Sequence ID")
  {
    node_response$Tables[[1]]$ColumnDescriptions[[nextCol]] <- node_args$Tables[[1]]$ColumnDescriptions[[i]]
    nextCol <- nextCol + 1
  }else if(grepl("Abundances Normalized", node_args$Tables[[1]]$ColumnDescriptions[[i]]$ColumnName)){
    node_response$Tables[[1]]$ColumnDescriptions[[nextCol]] <- node_args$Tables[[1]]$ColumnDescriptions[[i]]
    
    node_response$Tables[[1]]$ColumnDescriptions[[nextCol]]$ColumnName <- names(dat_filtered)[nextCol]
    node_response$Tables[[1]]$ColumnDescriptions[[nextCol]]$Options[[1]] <- "AbundancesFiltered"
    nextCol <- nextCol + 1
  }
}


toJSON(node_response, indent = 2) |>
  cat(file = node_args$ExpectedResponsePath)
