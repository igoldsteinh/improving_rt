# improving_rt
This repository has all code and data necessary to generate the results reported in Incorporating Testing Volume into  Estimation of Effective Reproduction Number Dynamics. All code was written in R or Stan. Each file that produces visualizations has the names of all files used to create that visualization. 

## Navigation

### Simulated Results
Files used to generate simulated outbreak data are located in the [simulation scripts folder](https://github.com/igoldsteinh/improving_rt/tree/main/R/simulation_studies/simulation_scripts). Files used to analyze simulated outbreak data are located in the [simulation analysis scripts folder](https://github.com/igoldsteinh/improving_rt/tree/main/R/simulation_studies/simulation_analysis_scripts). Much of this code was written to run on a computing cluster using the [Slurm scheduler](https://slurm.schedmd.com/slurm.html).  Pre-processed results and code for creating the manuscript figures are stored in the [simulation results folder](https://github.com/igoldsteinh/improving_rt/tree/main/R/simulation_studies/simulation_results). 


### Real Data Analysis
Data from the [California Open Data Portal](https://data.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state) used in this study is stored in the [data folder](https://github.com/igoldsteinh/improving_rt/tree/main/data). Code for analyzing this data is available in the [real data analysis scripts folder](https://github.com/igoldsteinh/improving_rt/tree/main/R/realdata_analysis/realdata_analysis_scripts). Again, much of this code was written to be run on a computing cluster. Pre-processed results and code for creating the manuscript figures are stored in the [real data results folder](https://github.com/igoldsteinh/improving_rt/tree/main/R/realdata_analysis/realdata_results). 