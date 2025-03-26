# test-sample_data.R
# test a few sample data sets


# filter protein abundance statistics for a file with a single replicate
test_that("Single replicate works", {

  #### Read in node args and data from PD ####

  # read in node_args
  node_args <- read_sample_node_args(0)

  expect_type(node_args, 'list')
  expect_length(node_args$Tables, 1)


  # fetch column types
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


  #### Write node response file ####

  node_response <- write_PD_response(node_args, dat_filtered, colTypes,
                                     startsWith = c('Abundance.Ratio', 'Abundances.Grouped'),
                                     save = FALSE)

  expect_type(node_response, 'list')
  expect_length(node_response$Tables, 1)
  expect_equal(node_response$Tables[[1]]$TableName, 'Proteins')

  expect_equal(node_response$Tables[[1]]$ColumnDescriptions[[2]]$ColumnName,
               "Filtered.Abundances.Grouped_F7")

  expect_equal(node_response$Tables[[1]]$ColumnDescriptions[[4]]$ColumnName,
               "Filtered.Abundance.Ratio_F7/F8")
})


# filter protein abundance statistics for a file with multiple replicates
test_that("Multiple replicates works", {

  #### Read in node args and data from PD ####

  # read in node_args
  node_args <- read_sample_node_args(1)

  expect_type(node_args, 'list')
  expect_length(node_args$Tables, 1)


  # fetch column types
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


  #### Write node response file ####

  node_response <- write_PD_response(node_args, dat_filtered, colTypes,
                                     startsWith = c('Abundance.Ratio', 'Abundances.Grouped'),
                                     save = FALSE)

  expect_type(node_response, 'list')
  expect_length(node_response$Tables, 1)
  expect_equal(node_response$Tables[[1]]$TableName, 'Proteins')

  expect_equal(node_response$Tables[[1]]$ColumnDescriptions[[2]]$ColumnName,
               "Filtered.Abundances.Grouped_Uninfected.Female")

  expect_equal(node_response$Tables[[1]]$ColumnDescriptions[[8]]$ColumnName,
               "Filtered.Abundance.Ratio_Uninfected.Female/Uninfected.Male")
})
