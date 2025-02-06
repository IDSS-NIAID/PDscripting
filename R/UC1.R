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

# transform template into valid JSON string 
# a) set CurrentWorkflowID
responseJSON <- sub('\\$CWFID\\$', PD_json_in$CurrentWorkflowID, NR_Template)

# b) set path of output tables
nodeResponsePath <- dirname(PD_json_in$ExpectedResponsePath)
responseJSON2 <- gsub("\\$PATH\\$", nodeResponsePath, responseJSON)
# print(responseJSON2)

# write node_response.json
jsonOutFile <- PD_json_in$ExpectedResponsePath
# print(jsonOutFile)
jsonfileconn <- file(jsonOutFile)
writeLines(responseJSON2, jsonfileconn)
close (jsonfileconn)

# --------------------------- PHASE 2 -------------------------------

# read in the node_response.json file (the one just written out)
NR_json_in <- fromJSON(file=jsonOutFile)
# print (NR_json_in)

library(data.table)

# read in input table
# print (PD_json_in$Tables[[1]]$DataFile)
protein.input <- data.table( 
  'Proteins Unique Sequence ID'=character(), 
  'Coverage in Percent'=numeric()
)
protein.input <- fread(PD_json_in$Tables[[1]]$DataFile, integer64 = "character", header = TRUE)
# print('Proteins data loaded')

#
#
#	instead of a non-R style 
#		(create new table then loop through the input table adding new table rows):
###
###	# create new table 
###	protein.output = data.table(
###	  'Proteins Unique Sequence ID'=character(),
###	  'Good Coverage'=logical()
###	)
###
###	# Loop through proteins
###	for (row in 1:nrow(protein.input)) {
###	  ins <- protein.input[row]
###	  ID = ins$'Proteins Unique Sequence ID'
###
###	  if (ins$'Coverage in Percent' > 20) {
###	     coverage = TRUE
###	  } else {
###	     coverage = FALSE
###	  }
###	  temprow = data.table(
###				'Proteins Unique Sequence ID'=ID,
###				'Good Coverage'=coverage
###				)
###	  protein.output <- rbind(protein.output, temprow)
###	}
###
### # Write output table
### # print (NR_json_in$Tables[[1]]$DataFile)
###	write.table(protein.output, file = NR_json_in$Tables[[1]]$DataFile, sep='\t', row.names=FALSE)
###
#	do the same in R style 
#		(use ifelse to create a new column from input table data, 
#		add it to the input table, 
#		remove no longer needed column from the input table,
#		and write thus modified input table as the output table):
#

# create a new logical column and add in to the input table
coverage <- ifelse(protein.input$'Coverage in Percent' > 20, TRUE, FALSE)
protein.input$'Good Coverage' <- coverage

#remove 2nd column from the input table
protein.input[,'Coverage in Percent':=NULL]


# Write modified input table as output table
# print (NR_json_in$Tables[[1]]$DataFile)
write.table(protein.input, file = NR_json_in$Tables[[1]]$DataFile, sep='\t', row.names=FALSE)

print('Done.')
