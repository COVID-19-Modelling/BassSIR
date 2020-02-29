# BassSIR
Infectious disease modelling with Bass SIR models


## Install package
In R, you can install from github
```r
devtools::install_github("TimeWz667/BassSIR", upgrade = FALSE)
```

Load the package before use it
```r
library(BassSIR)
```

## Contents

- [Modelling](#modelling)
  - [Data prepartion](#data-preparation)
  - [Model fitting](#model-fitting)
  - [Simulation](#simulation)
  - [Model comparison](#model-comparison)
  - [Scenario analysis](#scenario-analysis)
  
- [Output and visualisation](#output-and-visualisation)

- [Academic contact](#academic-contact)



## Modelling
### Data preparation
```r
cases <- as_bass_data(n_covid19$Hubei, id = "Hubei")
est <- BassSIR::fit(cases, r_rec = 1/22.2, r_death = 1/22.3, type = "BassSIR")
sim <- simulate(est, nsim = 1000)
```

### Model fitting
Fit BassSIR model to data
```r
est_bass <- BassSIR::fit(cases, r_rec = 1/22.2, r_death = 1/22.3, type = "BassSIR")
summary(est_bass)
```

Fit SIR model to data
```r
est_sir <- BassSIR::fit(cases, r_rec = 1/22.2, r_death = 1/22.3, type = "SIR")
summary(est_sir)
```

### Model comparison
```r
compare_models(BassSIR = est_bass, SIR = est_sir)
```

### Simulation
```r
sim <- simulate(est_bass, nsim = 1000)
```

### Scenario analysis
```r
zero_kappa <- function(pars) {
  pars$kappa <- rep(0, length(pars$kappa))
  return(pars)
}
lockdown <- run_scenario(sim, zero_kappa)
```

**compare_scenarios** will output two group of time-series

- **Trajectories** The trends of simulations of each variable
- **Change** Calculate the value changes from the baseline case 

```{r}
cp <- compare_scenarios(sim, Lockdown = lockdown, fn_change = "Yt")
```

## Output and visualisation
```r
TBA
```


## Academic contact

Chu-Chang Ku,

Health Economic and Decision Science, University of Sheffield, UK

Email: c.ku@sheffield.ac.uk


## License

MIT (c) 2020 Chu-Chang Ku, Ta-Chou Ng, and Hsien-Ho Lin

