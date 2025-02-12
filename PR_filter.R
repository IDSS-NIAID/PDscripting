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

if(!interactive())
{
  # this will run on the scripting node (it runs in batch mode)
  node_args <- fromJSON(file = commandArgs(trailingOnly = TRUE)[1])
}else{
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
                   new_name = case_when(grepl("Abundances Normalized", name) ~ str_replace(name, ' ', '.') |> # convert to "Abundances.Normalized <sample Id>.Sample"
                                                                               str_replace(' ', '_') |>       # convert to "Abundances.Normalized_<sample Id>.Sample"
                                                                               str_replace_all(' ', '.'),     # convert to "Abundances.Normalized_<sample.Id>.Sample"

                                         grepl("Found in Sample", name)      ~ str_replace(name, ' ', '.') |> # convert to "Found.in.Samples <sample Id>.Sample"
                                                                               str_replace(' ', '.') |> 
                                                                               str_replace(' ', '_') |>       # convert to "Found.in.Samples_<sample Id>.Sample"
                                                                               str_replace_all(' ', '.'),     # convert to "Found.in.Samples_<sample.Id>.Sample"

                                        grepl("Abundance Ratio", name)       ~ str_replace(name, ' ', '.') |> # convert to "Abundance.Ratio* <sample Id>  <sample Id>"
                                                                               str_replace('Ratio log2', 'Ratio.log2') |>
                                                                               str_replace(' ', '_') |>       # convert to "Abundance.Ratio*_<sample Id>  <sample Id>"
                                                                               str_replace('  ', '/') |>      # convert to "Abundance.Ratio*_<sample Id>/<sample Id>"
                                                                               str_replace(' ', '.'),         # convert to "Abundance.Ratio*_<sample.Id>/<sample.Id>"
                                        
                                        grepl('Abundances Grouped', name)    ~ str_replace(name, ' ', '.') |> # convert to "Abundances.Grouped <sample Id>"
                                                                               str_replace(' ', '_') |>       # convert to "Abundances.Grouped_<sample Id>"
                                                                               str_replace_all(' ', '.'),     # convert to "Abundances.Grouped_<sample.Id>"
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


#### Filtering checks ####

# convert to long format for summarization
dat_check <- dat |>
  select(`Proteins Unique Sequence ID`,
         Accession, 
         starts_with("Abundances.Normalized"),
         starts_with("Found.in.Sample")) |>
  
  pivot_longer(cols = starts_with("Abundances.Normalized") | starts_with("Found.in.Sample"),
               names_to = c(".value", "sample_replicate"),
               names_pattern = "(.*)_(.*)") |>
  
         # pull out the condition
  mutate(condition = str_replace(sample_replicate, 'F\\d+\\.Sample', '') |>
                     str_replace('BR\\d+\\.', '') |>
                     str_replace('^\\.', ''),
         
         # if there is no condition to group by, simply go with the sample replicate
         sample = ifelse(condition == '', 
                         str_replace(sample_replicate, '\\.Sample', ''), 
                         condition)) |>

  # filter out samples with too many missing or low values
  # count good abundances per sample
  group_by(`Proteins Unique Sequence ID`, Accession, sample) |>
  mutate(thresh_good = ceiling(n() * prop_good),
         ngood = sum(!is.na(Abundances.Normalized) & Abundances.Normalized > llod)) |>
  ungroup() |>
  
  # filter any samples that don't meet the threshold
  mutate(drop = ngood >= thresh_good,
         Abundances.Normalized = ifelse(drop, Abundances.Normalized, NA)) |>
  
  # drop unnecessary columns
  select(-Accession, -ngood, -thresh_good, -sample_replicate, -condition) |>

  # pivot back to wide format
  pivot_wider(names_from = sample,
              values_from = c("drop", "Abundances.Normalized", "Found.in.Sample"))
  

#### Do the actual filtering ####

# convert ratio columns to long format
dat_filtered <- dat |>
  
  select(`Proteins Unique Sequence ID`,
         starts_with("Abundance.Ratio")) |>
  
  pivot_longer(cols = starts_with("Abundance.Ratio"),
               names_to = c(".value", "samples"),
               names_pattern = "(.*)_(.*)") |>
  
  mutate(numer = paste0('drop_', str_split(samples, '\\/', simplify = TRUE)[, 1]),
         denom = paste0('drop_', str_split(samples, '\\/', simplify = TRUE)[, 2]),
         drop = FALSE)

# check against filtering checks
for(i in 1:nrow(dat_filtered))
{
  # row of dat_check to check against
  j <- which(dat_check$`Proteins Unique Sequence ID` == dat_filtered$`Proteins Unique Sequence ID`[i])
  
  # check if numerator and denominator are both good
  if(!dat_check[[ dat_filtered$numer[i] ]][j] |
     !dat_check[[ dat_filtered$denom[i] ]][j])
  {
    dat_filtered$drop[i] <- TRUE
  }
}

# filter and make wide again
dat_filtered <- dat_filtered |>
  filter(!drop) |>
  select(-drop, -numer, -denom) |>
  
  pivot_wider(names_from = samples,
              values_from = starts_with("Abundance.Ratio"))

# update column names
colnames(dat_filtered) <- str_replace(colnames(dat_filtered), 'Abundance', 'Filtered.Abundance')


#### Write output file ####

out_dir <- dirname(node_args$Tables[[1]]$DataFile)
write_delim(dat_filtered, file = file.path(out_dir, 'to_pd.tsv'), delim = '\t')


#### Write node response file ####

node_response = list(CurrentWorkflowID = node_args$CurrentWorkflowID,
                     Tables = list(list(TableName = "Proteins",
                                        DataFile  = file.path(out_dir, 'to_pd.tsv'),
                                        DataFormat = "CSV",
                                        Options = NULL,
                                        ColumnDescriptions = list())
                                  )
                    )

# update colTypes with new names
colTypes_updt <- mutate(colTypes,
                        new_name = str_replace(new_name, fixed('Abundance.Ratio'), 'Filtered.Abundance.Ratio'))


# loop over new column names
nextCol <- 1
for(nn in names(dat_filtered))
{
  # figure out which old name we are looking for
  on <- colTypes_updt$name[colTypes_updt$new_name == nn]
  
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


toJSON(node_response, indent = 2) |>
  str_replace_all(fixed('"{}"'), '{}') |> # for some reason it wants to quote empty lists
  str_replace_all(fixed("[\n\n          ]\n"), "{}") |> # for some reason we end up with these, too
  cat(file = node_args$ExpectedResponsePath)
