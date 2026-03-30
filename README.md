########################################################################################################################################################################
# Causal Inference: Simulation and real data application implementation code instructions
########################################################################################################################################################################
## How to Run

### Simulation (ATE)
- Go to: `Final code/Simulation/`
- Run: `ATE_simulation_run.R`

### Real Data Analysis
- Go to: `Final code/Real_data_code/`
- Run: `ATE_real_data_analysis_final code.R`

### Missing Covariates Simulation
- Go to: `Covariates_missing/Final simulation code/`
- Run: `Covariates_missing_final_code.R`

---

## Notes

- Helper files are automatically used by the main scripts  
- No need to run helper files separately  
- Ensure required R packages are installed  





This repository provides code to reproduce both simulation studies and real data analyses for Average Treatment Effect (ATE) estimation in causal inference. 
The main scripts for each analysis are located in the corresponding folders, while additional R files provide supporting helper functions that are automatically used by the main scripts. Below is the final flow chart:
ATE_estimation/
│
├── Final code/
│   │
│   ├── Simulation/
│   │   ├── ATE_simulation_run.R                 # Main file to generate simulation results and boxplots
│   │   └── ate_functions.R                      # Helper functions for simulation
│   │
│   └── Real_data_code/
│       ├── ATE_real_data_analysis_final code.R  # Main file for real data analysis
│       └── ATE_real_data_analysis_functions.R   # Helper functions for real data




########################################################################################################################################################################
# Covariates Missing: Simulation implementation code instructions
########################################################################################################################################################################
This repository provides code to reproduce simulation studies table missing covariates. 
The main scripts for each analysis are located in the corresponding folders, while additional R files provide supporting helper functions that are automatically used by the main scripts. Below is the final flow chart:
Covariates_missing/
│
├── Final simulation code/
│   │
│   │   ├── Covariates_missing_final_code.R      # Main file to generate simulation results and table
│   │   └── calibration_with_missing_covariates_functions.R     # Helper functions for simulation
