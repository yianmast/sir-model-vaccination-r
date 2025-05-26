
# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

# MODEL INPUTS:

# Vector storing the initial number in each compartment (at timestep 0)

N <- 300000

initial_state_values <- c(S = 0.5*N,  
                          I = 0.05*N,      
                          R = 0.45*N)     
# the exact proportions here don't matter, we have chosen an infection prevalence of 5% to start with in line with the
# information that the disease is thought to be relatively rare in this population. 
# As described above, any initial conditions will converge on the same endemic equilibrium, given enough time.

# Vector storing the parameters describing the transition rates in units of years^-1
parameters <- c(beta = 365/1,            # the infection rate
                gamma = 365/20,          # the rate of recovery
                mu = 1/3,                # the background mortality rate
                b = 1/3,                 # the birth rate
                p_vacc = 0,              # the neonatal vaccine coverage
                sigma = 0)               # the rate at which immunity wanes                   

# TIMESTEPS:

# Vector storing the sequence of timesteps to solve the model at
times <- seq(from = 0, to = 5, by = 1/365)   # from 0 to 5 years in daily intervals
# We are simulating over a period of 10 years to allow the model to come to equilibrium. 
# You might need a different timespan depending on the initial conditions you chose.

# SIR MODEL FUNCTION: 

# The model function takes as input arguments (in the following order): time, state and parameters
sir_model <- function(time, state, parameters) {  

    with(as.list(c(state, parameters)), {  # tell R to unpack variable names from the state and parameters inputs   
        
    # Calculating the total population size N (the sum of the number of people in each compartment)
      N <- S+I+R
      
    # Defining lambda as a function of beta and I:
      lambda <- beta * I/N
        
    # The differential equations
      dS <- -lambda*S - mu*S + (1-p_vacc)*b*N + sigma*R           
      dI <- lambda*S - gamma*I  - mu*I           
      dR <- gamma*I - mu*R  + p_vacc*b*N - sigma*R   
    
    # Because this is a neonatal vaccine (given straight after birth), we model this simply as a proportion p_vacc
    # of births entering the R compartment (Recovered/Immune),
    # with the remaining births (1-p_vacc) entering the susceptible compartment.
      
    # Output the number in the S, I and R compartments at each timestep 
    # (in the same order as the input state variables)
    return(list(c(dS, dI, dR))) 
    })
  
}

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))

# PLOT

output_long <- melt(as.data.frame(output), id = "time")                  # turn output dataset into long format

# Calculating the proportion in each compartment as a column in the long-format output
output_long$proportion <- output_long$value/sum(initial_state_values)

# Plot the proportion in the S, I and R compartments over time
ggplot(data = output_long,                                               # specify object containing data to plot
       aes(x = time, y = proportion, colour = variable, group = variable)) +  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Baseline prevalence of susceptible, infected and recovered animals over time")

# The prevalence seems to have stabilised by 2 years:
# calculating the prevalence in year 2
output_long$proportion[round(output_long$time,0) == 2 & output_long$variable == "I"][1]

# Note that here we are selecting the proportion infected at timestep 2, 
# but since the timesteps are not exact numbers, we are selecting the timesteps that, when rounded to 0 decimals,
# is and display only the first one of those

output[output$time == 2,]  # print state values at equilibrium

# Before introducing a vaccine, change initial state values to baseline endemic equilibrium
initial_state_values <- c(S = 15254,  
                          I = 5125,      
                          R = 279621)    

parameters["p_vacc"] <- 0.5   # try different coverage values to find an endemic prevalence of half the baseline prevalence

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                             times = times, 
                             func = sir_model,
                             parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")                  # turn output dataset into long format

# Calculating the proportion in each compartment as a column in the long-format output
output_long$proportion <- output_long$value/sum(initial_state_values)

# PLOT

# Plot the proportion in the I compartment over time: divide the number in I over time by the total population size
ggplot(data = output,                                                    # specify object containing data to plot
       aes(x = time, y = output$I/N)) +                                  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Infection prevalence with neonatal vaccine coverage of 50%")

# Calculating the prevalence in year 5 (introduction of the vaccine at first perturbs the initial equilibrium, but
# we are interested in the new endemic equilibrium achieved with vaccination)
output_long$proportion[round(output_long$time,0) == 5 & output_long$variable == "I"][1]

parameters["p_vacc"] <- 0.95   # try different coverage values to see if the disease persists

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                             times = times, 
                             func = sir_model,
                             parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")                  # turn output dataset into long format

# Calculating the proportion in each compartment as a column in the long-format output
output_long$proportion <- output_long$value/sum(initial_state_values)

# PLOT

# Plot the proportion in the I compartment over time: divide number in I over time by the total population size
ggplot(data = output,                                                    # specify object containing data to plot
       aes(x = time, y = output$I/N)) +                                  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Infection prevalence with neonatal vaccine coverage of 95%")

# Double-check how many animals remain infected at the 5 year timestep
output[round(output$time,0) == 5,"I"][1]

# BASELINE SCENARIO WITH WANING IMMUNITY

initial_state_values <- c(S = 0.5*N,  
                          I = 0.05*N,      
                          R = 0.45*N)

parameters["p_vacc"] <- 0
parameters["sigma"] <- 1

waning_baseline <- as.data.frame(ode(y = initial_state_values, 
                                     times = times, 
                                     func = sir_model,
                                     parms = parameters))

waning_baseline_long <- melt(as.data.frame(waning_baseline), id = "time")   

# Calculating the proportion in each compartment
waning_baseline_long$proportion <- waning_baseline_long$value/sum(initial_state_values)

# VACCINE SCENARIO WITH WANING IMMUNITY

parameters["p_vacc"] <- 0.5
parameters["sigma"] <- 1

waning_vacc <- as.data.frame(ode(y = initial_state_values, 
                                     times = times, 
                                     func = sir_model,
                                     parms = parameters))

waning_vacc_long <- melt(as.data.frame(waning_vacc), id = "time")   

# Calculating the proportion in each compartment
waning_vacc_long$proportion <- waning_vacc_long$value/sum(initial_state_values)

# PLOTTING THE 2 SCENARIOS

# Baseline
ggplot(data = waning_baseline,                                           # specify object containing data to plot
       aes(x = time, y = I/N)) +                                         # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Baseline prevalence with waning immunity") +
  ylim(c(0,0.4))                                                         # define y axis range to be plotted (0-0.4)

# With vaccine
ggplot(data = waning_vacc,                                               # specify object containing data to plot
       aes(x = time, y = I/N)) +                                         # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Prevalence with waning immunity and vaccine coverage of 50%") +
  ylim(c(0,0.4))                                                         # define y axis range to be plotted (0-0.4)

# Calculating the baseline prevalence with waning immunity
waning_baseline_prev <- waning_baseline_long$proportion[round(waning_baseline_long$time,0) == 2 & 
                                                        waning_baseline_long$variable == "I"][1]

# Calculating the endemic prevalence with waning immunity and neonatal vaccination coverage of 50%
waning_vacc_prev <- waning_vacc_long$proportion[round(waning_vacc_long$time,0) == 2 & 
                                                waning_vacc_long$variable == "I"][1]

# Calculating the reduction in prevalence achieved with 50% neonatal vaccine coverage:
1-waning_vacc_prev/waning_baseline_prev

# VACCINE SCENARIO WITH WANING IMMUNITY: INCREASING COVERAGE TO 100%

parameters["p_vacc"] <- 1
parameters["sigma"] <- 1

waning_vacc <- as.data.frame(ode(y = initial_state_values, 
                                     times = times, 
                                     func = sir_model,
                                     parms = parameters))

# With vaccine
ggplot(data = waning_vacc,                                      # specify object containing data to plot
       aes(x = time, y = I/N)) +  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (years)")+                                                  # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title  
       title = "Prevalence with waning immunity and vaccine coverage of 100%") +
ylim(c(0,0.4))

# BASELINE SCENARIO WITH WANING IMMUNITY

parameters["p_vacc"] <- 0
parameters["sigma"] <- 1/2.5

waning_baseline <- as.data.frame(ode(y = initial_state_values, 
                                     times = times, 
                                     func = sir_model,
                                     parms = parameters))

waning_baseline_long <- melt(as.data.frame(waning_baseline), id = "time")   

# Calculating the proportion in each compartment
waning_baseline_long$proportion <- waning_baseline_long$value/sum(initial_state_values)

# VACCINE SCENARIO WITH WANING IMMUNITY

parameters["p_vacc"] <- 1
parameters["sigma"] <- 1/2.5

waning_vacc <- as.data.frame(ode(y = initial_state_values, 
                                     times = times, 
                                     func = sir_model,
                                     parms = parameters))

waning_vacc_long <- melt(as.data.frame(waning_vacc), id = "time")   

# Calculating the proportion in each compartment
waning_vacc_long$proportion <- waning_vacc_long$value/sum(initial_state_values)

# Calculating the baseline prevalence with slower waning immunity
print("Baseline prevalence:")
waning_baseline_long$proportion[round(waning_baseline_long$time,0) == 2 & waning_baseline_long$variable == "I"][1]
# Calculating the endemic prevalence with slower waning immunity and neonatal vaccination coverage of 100%
print("Prevalence with neonatal vaccine coverage of 100%:")
waning_vacc_long$proportion[round(waning_vacc_long$time,0) == 2 & waning_vacc_long$variable == "I"][1]










