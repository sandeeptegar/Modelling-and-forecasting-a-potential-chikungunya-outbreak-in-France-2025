# Instructions to Run the Code

## 1. *Ae. albopictus* Population Dynamics

**Main code file:** albopictus_pop.ipynb  
**Supporting function file:** AlbopictusFun.jl  

### Output:  
This script reproduces the population dynamics of adult female *Aedes albopictus* and their oviposition activity.  

### Required Data:  
- Data files: `AMgam.csv`, `Car_Tol_1.csv`, `Fin_LSurv.csv`, `LDgam.csv`, `WL_re.csv`  
- Climate file: `latitude_longitude.csv`  
  (For Fano: latitude = 43.82, longitude = 13.026)
- Source: [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview)   

---

## 2. Chikungunya Outbreak Simulation

**Main code file:** bernis_base_2025.ipynb  
**Supporting function files:** ChikungunyaFunSR1.jl, seeding_run.jl, read_case_files.jl 

### Required Data:  
- Data files: `AMgam.csv`, `Car_Tol_1.csv`, `Fin_LSurv.csv`, `LDgam.csv`, `WL_re.csv`  
- Climate file: `latitude_longitude.csv`  
  (For Bernis: latitude = 43.77, longitude = 4.29) 

---

## 3. Climate scenarios

**Main code file:** bernis_clim_filled.ipynb  

### Output:  
This script calculates the ten climate scenarios based on the historic climate for 2015–2024.  

### Required Data:  
- Climate file: `Bernis_climate_1979_2025.csv` (Observed climate data for Bernis up to the latest date)
- Source: [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview)     

---

## 4. How to Reproduce Plots for Each Scenarios

### Scenario 2015:  
- Use the created scenario file for the year 2015 and write in the format: `latitude_longitude.csv`  
- use this climate file to run the code
- repeat the same for each climate scenario from 2015-2024 

### Final plot  
- After running the code for each climate scenario, find mean, 95% CI, min, and max from the 10 output files obtained from 10 climate scenarios.  
- Finally plot the above information, i.e., mean, mean, 95% CI, min, and max.
- The method applies for both mosquito population dynamics and transmission dynamics.

---

### Required Data:  
- Population data per sq. km around Bernis city center  
- Source: [Eurostat - Population density data](https://ec.europa.eu/eurostat/statistics-explained/index.php?oldid=596753)  
