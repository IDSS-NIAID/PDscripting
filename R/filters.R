# filters.R
# functions for filtering data

#' PD_check_missingness
#' 
#' Check for missingness in data by condition and flag conditions with too many missing or low values.
#' 
#' @param dat data.frame to check
#' @param startsWith Character value identifying columns to check (passed to `starts_with()`)
#' @param prop_good Numeric proportion of good values to keep
#' @param llod Numeric lower limit of detection
#' 
#' @return data frame with missingness information and which rows to drop
#' @export
#' @importFrom dplyr filter group_by mutate select ungroup
#' @importFrom rlang sym
#' @importFrom stringr str_extract str_replace
#' @importFrom tidyr pivot_longer pivot_wider
PD_check_missingness <- function(dat, startsWith, prop_good = 0.5, llod = 1e6)
{
  # convert to long format for summarization
  dat_check <- dat |>
    select(`Proteins Unique Sequence ID`,
           Accession, 
           starts_with(startsWith)) |>
    
    pivot_longer(cols = -(`Proteins Unique Sequence ID`:Accession),
                 names_to = c(".value", "sample_replicate"),
                 names_pattern = "(.*)_(.*)") |>
    
    # pull out the condition
    mutate(condition = str_replace(sample_replicate, 'F\\d+\\.Sample', '') |>
             str_replace('BR\\d+\\.', '') |>
             str_replace('^\\.', ''),
           
           # if there is no condition to group by, simply go with the sample replicate
           sample = ifelse(condition == '', 
                           str_replace(sample_replicate, '\\.Sample', ''), 
                           condition),
           
           # get the replicate
           rep = str_extract(sample_replicate, 'F\\d+')) |>
    
    # filter out samples with too many missing or low values
    # count good abundances per sample
    group_by(`Proteins Unique Sequence ID`, Accession, sample) |>
    mutate(thresh_good = ceiling(n() * prop_good),
           ngood = sum(!is.na(!!sym(startsWith)) & !!sym(startsWith) > llod)) |>
    ungroup() |>
    
    # filter any samples that don't meet the threshold
    mutate(drop = ngood < thresh_good,
           Abundances.Normalized = ifelse(drop, NA, Abundances.Normalized),
           id = paste(rep, sample, sep = '_')) |>
    
    # drop unnecessary columns
    select(-Accession, -ngood, -thresh_good, -sample_replicate, -condition, -starts_with("Found.in.Sample"),
           -rep, -sample) |>
    
    # pivot back to wide format
    pivot_wider(names_from = id,
                values_from = c("drop", "Abundances.Normalized"))
  
  colnames(dat_check) <- str_replace(colnames(dat_check), '_(.*?)_', '_') # fix annoying duplication of labels in some cases
  
  return(dat_check)
}


#' PD_filter
#' 
#' Filter data by condition
#' 
#' @param dat data.frame to filter
#' @param dat_check data.frame with missingness information from `PD_check_missingness()`
#' @param startsWith Character value identifying columns to check (passed to `starts_with()`)
#' @param ratio Logical value indicating whether the columns to be filtered are ratios
#' 
#' @return data frame with filtered data
#' @export
PD_filter <- function(dat, dat_check, startsWith, ratio)
{
  # convert dat columns to long format
  dat_filtered <- dat |>
    
    select(`Proteins Unique Sequence ID`,
           starts_with(startsWith)) |>
    
    pivot_longer(cols = starts_with(startsWith),
                 names_to = c(".value", "samples"),
                 names_pattern = "(.*)_(.*)")
  
  if(ratio)
  {
    # if filtering by ratio, we need to check both numerator and denominator
    dat_filtered <- dat_filtered |>
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
  }else{
    # if filtering by grouped abundance, we only need to check the grouped abundance
    dat_filtered <- dat_filtered |>
      mutate(smp = paste0('drop_', samples),
             drop = FALSE)
    
    # check against filtering checks
    for(i in 1:nrow(dat_filtered))
    {
      # row of dat_check to check against
      j <- which(dat_check$`Proteins Unique Sequence ID` == dat_filtered$`Proteins Unique Sequence ID`[i])
      
      # check if sample is good
      if(!dat_check[[ dat_filtered$smp[i] ]][j])
      {
        dat_filtered$drop[i] <- TRUE
      }
    }
  }
  
  # filter and make wide again
  dat_filtered <- dat_filtered |>
    filter(!drop) |>
    select(`Proteins Unique Sequence ID`, samples, starts_with(startsWith)) |>
    
    pivot_wider(names_from = samples,
                values_from = starts_with(startsWith),
                names_glue = "{.value}_{.name}") # make sure we get the full names - this will create some annoying duplicate labels, but we'll fix that below
  
  # update column names
  colnames(dat_filtered) <- str_replace(colnames(dat_filtered), startsWith, paste0('Filtered.', startsWith)) |>
    str_replace("_(.*?)_", "_") # fix annoying duplication of labels in some cases
  
  return(dat_filtered)
}
