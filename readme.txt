This file contains notes on the contents and instructions to run the regionalmalariaweathermodels project

Run the whole regionalweatherapproaches.R script to generate simulation output images. This process takes **approximately 45 minutes** to generate images, the duration of each stage in the process are provided at the end of the simulation for reference.

The script contains packages necessary to pre-processing, simulation, post-processing, and themes for plot visualization.

# preprocessing

Initial conditions and the model parameter values for the vector and host populations are read from the vector_host_states_parameters.xlsx file.

Monthly temperature and rainfal data are read from era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era_x0.25_kano/benisha/limpopo.xlsx files then interpolated into daily values.

# simulation

The model function vectorhostModel contains the parameter functionals for Approach A, B and C as well as model equations.

The sim_model function simulates the model for different intervention coverages for each approach.

Results from parallel simulations of model scenarios (different vector control intervention coverage) for each approach are combined into the *result* object which contains all the simulation results for the 1950 to 2022 period.

# postprocessing

The convert_to_longer function splits model output results for each approach by patches (climatic regions) the total populations of humans and mosquitoes together with the cummulative incidence are calculated. These are then bound into the *output* object. 

Annual cummulative incidence predicted for each approach are determined using the incidence_summary functions.

# visualization

The following plots are generated from model outputs:

1. Mean monthly temperature/rainfall trends 
2. Aquatic (eggs, larvae, and pupae) mosquito populations in all approaches (columns) across climate regions (rows)
3. Adult (susceptible, exposed, and infectious) mosquito populations in all approaches (columns) across climate regions (rows)
4. Human (susceptible, exposed, asymptomatic, mild symptomatic, severe symptomatic, treated mild cases, treated severe cases, and recovered) populations predicted in all approaches (columns) across climate regions (rows)
5/6. Daily/Annual malaria incidence predicted in all approaches (columns) across climate regions (rows)

# comments

The number of cores used for parallel simulation and durations of installing packages, preprocessing, simulation, postprocessing, plotting, and the total time taken to run the whole script are provided at the end of the simulation.








