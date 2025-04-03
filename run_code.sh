#!/bin/bash

# Set the path to the R script
R_SCRIPT_PATH="./code/code.R"

# Run the R script
Rscript $R_SCRIPT_PATH

# Check if the R script ran successfully
if [ $? -eq 0 ]; then
    echo "R script completed successfully."
else
    echo "R script did not complete successfully."
fi


