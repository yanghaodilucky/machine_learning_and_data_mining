# Comparative Analysis of Clinical and Genomic Data for Predicting Breast Cancer Recurrence: Evaluating Machine Learning Models

## Description
The study utilized different combinations of datasets including clinical data, mRNA data, miRNA data, protein expression data, and genetic mutation data. Further, the study compared the effectiveness of XGBoost, Elastic Net, Decision Tree, Random Forest, and LASSO using clinical plus genetic mutation data.

## Getting Started
-After unzip this folder, using linux to get access to the project directory
cd MLDM_assessment_part2_2398736
ls
-Please make sure there are data folder, code folder, visullization folder and result folder. That is important because my result and plots will be generated there. And also renv folder to record the environment, a git_history.txt to see the past versons. 

-Then run the bash file: run_code.sh in MLDM_assessment_part2_2398736 folder
sh run_code.sh

## Data Visualization Approach
after running the code, there will be some plots in visulization folder and the tables in csv format in result format. 

## Results
By comparing the predictive abilities of different datasets, I found genetic mutation data significantly enhanced prediction accuracy. The results demonstrate that the XGBoost model performed optimally in handling complex genetic relationships, providing the highest prediction accuracy. 

### note: 
to make the folder small enough to submit, I make the git smaller and have a git_history.txt in git folder to demonstrate that I used git at the very first. 
