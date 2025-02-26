# sample_data.R
# This is a sample data file for testing purposes and examples.

library(dplyr)
library(readr)

root <- here::here()


##### Example 0 #####

# simple case where we only have one replicate per group
n <- 10

set.seed(92384)
dat <- tibble(`Proteins Unique Sequence ID` = round(runif(n) * -1e10, 0) |> as.character(),
              Accession = 'P51681',
              `Abundances Normalized F7 Sample` = runif(n) * 1e7,
              `Abundances Normalized F8 Sample` = runif(n) * 1e7,
              `Found in Sample F7 Sample` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F8 Sample` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Abundance Ratio F7  F8` = `Abundances Normalized F7 Sample` / `Abundances Normalized F8 Sample`,
              `Abundance Ratio log2 F7  F8` = log2(`Abundance Ratio F7  F8`),
              `Abundances Grouped F7` = `Abundances Normalized F7 Sample`,
              `Abundances Grouped F8` = `Abundances Normalized F7 Sample`)

# add some random spaces that PD seems to include in output
for(i in c(4,6,10))
{
  dat[[i]] <- paste0(' ', dat[[i]])
}

# add some missing values
dat$`Abundances Normalized F7 Sample`[1:3] <- NA
dat$`Abundances Normalized F8 Sample`[3:5] <- NA

write_delim(dat, file = file.path(root, 'inst', 'extdata', 'TargetProtein0.txt'),
            na = '', delim = '\t')


##### Example 1 #####

# more complicated case where we have 3 replicates per group
n <- 10

set.seed(92384)
dat <- tibble(`Proteins Unique Sequence ID` = round(runif(n) * -1e10, 0) |> as.character(),
              `Accession`  = 'P51681',
              `Abundances Normalized F4 Sample BR1 Uninfected Female` = runif(n) * 1e7,
              `Abundances Normalized F5 Sample BR2 Uninfected Female` = runif(n) * 1e7,
              `Abundances Normalized F6 Sample BR3 Uninfected Female` = runif(n) * 1e7,
              `Abundances Normalized F1 Sample BR1 Uninfected Male` = runif(n) * 1e7,
              `Abundances Normalized F2 Sample BR2 Uninfected Male` = runif(n) * 1e7,
              `Abundances Normalized F3 Sample BR3 Uninfected Male` = runif(n) * 1e7,
              `Found in Sample F4 Sample BR1 Uninfected Female` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F5 Sample BR2 Uninfected Female` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F6 Sample BR3 Uninfected Female` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F1 Sample BR1 Uninfected Male` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F2 Sample BR2 Uninfected Male` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Found in Sample F3 Sample BR3 Uninfected Male` = sample(c('Peak Found', 'High'), n, replace = TRUE),
              `Abundance Ratio Adj P-Value Uninfected Female  Uninfected Male` = runif(n),
              `Abundance Ratio P-Value Uninfected Female  Uninfected Male` = runif(n),
              `Abundance Ratio Uninfected Female  Uninfected Male` = mean(c(`Abundances Normalized F4 Sample BR1 Uninfected Female`,
                                                                             `Abundances Normalized F5 Sample BR2 Uninfected Female`,
                                                                             `Abundances Normalized F6 Sample BR3 Uninfected Female`)) /
                                                                      mean(c(`Abundances Normalized F1 Sample BR1 Uninfected Male`,
                                                                             `Abundances Normalized F2 Sample BR2 Uninfected Male`,
                                                                             `Abundances Normalized F3 Sample BR3 Uninfected Male`)),
              `Abundance Ratio log2 Uninfected Female  Uninfected Male` = log2(`Abundance Ratio Uninfected Female  Uninfected Male`),
              `Abundances Grouped Uninfected Female` = mean(c(`Abundances Normalized F4 Sample BR1 Uninfected Female`,
                                                              `Abundances Normalized F5 Sample BR2 Uninfected Female`,
                                                              `Abundances Normalized F6 Sample BR3 Uninfected Female`)),
              `Abundances Grouped Uninfected Male` = mean(c(`Abundances Normalized F1 Sample BR1 Uninfected Male`,
                                                            `Abundances Normalized F2 Sample BR2 Uninfected Male`,
                                                            `Abundances Normalized F3 Sample BR3 Uninfected Male`)),
              `Abundances Grouped Count Uninfected Female` = 3,
              `Abundances Grouped Count Uninfected Male` = 3)

# add some random spaces that PD seems to include in output
for(i in c(4:8, 10:14, 20, 22))
{
  dat[[i]] <- paste0(' ', dat[[i]])
}

# add some missing values
dat$`Abundances Normalized F4 Sample BR1 Uninfected Female`[1:3] <- NA
dat$`Abundances Normalized F5 Sample BR2 Uninfected Female`[2:4] <- NA
dat$`Abundances Normalized F6 Sample BR3 Uninfected Female`[3:5] <- NA

dat$`Abundances Normalized F1 Sample BR1 Uninfected Male`[3:5] <- NA
dat$`Abundances Normalized F2 Sample BR2 Uninfected Male`[4:6] <- NA
dat$`Abundances Normalized F3 Sample BR3 Uninfected Male`[5:7] <- NA

write_delim(dat, file = file.path(root, 'inst', 'extdata', 'TargetProtein1.txt'),
            na = '', delim = '\t')
