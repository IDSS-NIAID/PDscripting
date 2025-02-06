# R Script that reads in Proteins table and adds a new column to it.
#
# invoke as: r -f UC1.R --vanilla -q --args <full path of the node_args.json file>
# for example: r -f UC1.R --vanilla -q --args data/UC1/node_args.json

# The processing consists of two phases:
# 1. read node_args.json file and write node_response.json file
# 2. write output table


# Read command line arguments
arg <- commandArgs(trailingOnly = TRUE)
# print (arg)
# get the node_args.json file full path as the only arg after "--args"
inputFile <- arg[1]

if(FALSE)
{
  # when running manually, navigate to the R directory if we aren't already there
  root <- here::here()
  file.path(root, 'R') |>
    setwd()
  
  inputFile <- 'data/UC1/node_args.json'
}
# print (inputFile)

# --------------------------- PHASE 1 -------------------------------

# read node_args.json as JSON file
library(rjson)
PD_json_in <- fromJSON(file=inputFile)

# this is the node_response.json template that has
# two patterns ($CWFID$ and $PATH$) to be replaced 
# with the CurrentWorkflowID and the directory of 
# the ExpectedResponsePath correspondingly

NR_Template <- '
    {
      "CurrentWorkflowID": $CWFID$,
      "Tables": [
        {
          "TableName": "Proteins",
          "DataFile": "$PATH$/uc1_result.txt",
          "DataFormat": "CSV",
          "Options": {},
          "ColumnDescriptions": [
            {
              "ColumnName": "Proteins Unique Sequence ID",
              "ID": "ID",
              "DataType": "Long",
              "Options": {}
            },
            {
              "ColumnName": "Good Coverage",
              "ID": "",
              "DataType": "Boolean",
              "Options": {}
            }
          ]
        }
      ]
    }
'

sub('\\$CWFID\\$', PD_json_in$CurrentWorkflowID, NR_Template) |>
  gsub(pattern = "\\$PATH\\$",
       replacement = dirname(PD_json_in$ExpectedResponsePath)) |>
  cat(file = PD_json_in$ExpectedResponsePath)


# --------------------------- PHASE 2 -------------------------------

# read in the node_response.json file (the one just written out)
NR_json_out <- fromJSON(file=jsonOutFile)

library(dplyr)
library(readr)

# read in input table
protein.input <- read_delim(PD_json_in$Tables[[1]]$DataFile, delim = '\t') |>
  
  # create a new logical column and add in to the input table
  mutate(`Good Coverage` = ifelse(`Coverage in Percent` > 20, TRUE, FALSE)) |>

  #remove 2nd column from the input table
  select(-`Coverage in Percent`)


# Write modified input table as output table
write_delim(protein.input, file = NR_json_out$Tables[[1]]$DataFile, sep='\t', row.names=FALSE)
