# PR_fitler
# Filter replicates for valid protein ratios

# documentation on custom script integration: https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1325195659v1


#### Setup ####
#' Requested tables and columns: "Proteins: Accession, Abundances Normalized, Found in Samples"


library(rjson)
library(stringr)
library(PDscripting)

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
  root <- here::here()

  # for manual debugging
  node_args <- fromJSON(file = file.path(root, 'debug', 'node_args.json'))

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
colTypes <- PD_colTypes(node_args$Tables[[1]]$ColumnDescriptions)

# read in data
dat <- PD_clean(node_args$Tables[[1]]$DataFile, colTypes)


#### Filtering checks ####

dat_check <- PD_check_missingness(dat, startsWith = "Abundances.Normalized")


# actual filtering
grouped_filter <- PD_filter(dat, dat_check, "Abundances.Grouped", ratio = FALSE)
ratio_filter   <- PD_filter(dat, dat_check, "Abundance.Ratio",    ratio = TRUE)

# join filtered data
dat_filtered <- full_join(grouped_filter, ratio_filter,
                          by = "Proteins Unique Sequence ID")


#### Write output file ####

write_delim(dat_filtered,
            file = file.path(dirname(node_args$Tables[[1]]$DataFile),
                             'to_pd.tsv'),
            delim = '\t')


#### Write node response file ####

write_PD_response(node_args, dat_filtered, colTypes,
                  startsWith = c('Abundance.Ratio', 'Abundances.Grouped'))
