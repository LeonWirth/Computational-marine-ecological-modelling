library(deSolve)
library(ggplot2)
library(reshape2)

# Define the light function
light <- function(P ,D , Z, I0, kw, kp, dZ) {
  L <- cumsum(kp * (P+D) * dZ)  # Cumulative light attenuation due to phytoplankton
  I <- I0 * exp(-kw * Z - L)  # Light intensity at depth Z (micromol photons/m^2)
  return(I)
}

# Define the derivatives function
derivatives <- function(t, Y, parms) {
  # Split Y into N and P
  N <- Y[1:length(Z)]
  P <- Y[(length(Z) + 1):(2 * length(Z))]
  D <- Y[2*(length(Z))+1:(3 * length(Z))]
  
  # Calculate light intensity
  I <- light(P, D, Z, I0, kw, kp, dZ)
  
  # Define flux vectors
  JaP <- rep(0, length(Z) + 1)
  JdP <- rep(0, length(Z) + 1)
  
  JaN <- rep(0, length(Z) + 1)
  JdN <- rep(0, length(Z) + 1)
  
  JaD <- rep(0, length(Z) + 1)
  JdD <- rep(0, length(Z) + 1)
  
  # Define fluxes
  for (i in 2:length(Z)) {
    JaP[i] <- u * P[i - 1]
    JdP[i] <- -d * (P[i] - P[i - 1]) / dZ
    JaN[i] <- 0 * N[i - 1]
    JdN[i] <- -d * (N[i] - N[i - 1]) / dZ
    JaD[i] <- u * D[i - 1]
    JdD[i] <- -d * (D[i] - D[i - 1]) / dZ 
  }
  
  # Boundary conditions
  JaP[1] <- 0
  JdP[1] <- 0
  JaP[length(Z) + 1] <- 0
  JdP[length(Z) + 1] <- 0
  
  JaN[1] <- 0
  JdN[1] <- 0
  JaN[length(Z) + 1] <- 0
  JdN[length(Z) + 1] <- -d * (Nb - N[length(Z)]) / dZ
  
  JaD[1] <- 0
  JdD[1] <- 0
  JaD[length(Z) + 1] <- u * D[length(Z)] 
  JdD[length(Z) + 1] <- 0
  
  JP <- JaP + JdP
  JN <- JaN + JdN
  JD <- JaD + JdD
  
  dPdt <- rep(0, length(Z))
  dNdt <- rep(0, length(Z))
  dDdt <- rep(0, length(Z))
  
  for (i in 1:length(Z)) {
    growth_rate <- gmax * min(I[i] / (Hl + I[i]), N[i] / (Hn + N[i]))
    dNdt[i] <- (-(JN[i + 1] - JN[i]) / dZ) - a * growth_rate * P[i] + e * a * t * D[i]
    dPdt[i] <- (-(JP[i + 1] - JP[i]) / dZ) + (-m + growth_rate) * P[i]
    dDdt[i] <- (-(JD[i + 1] - JD[i]) / dZ) + m * P[i] - t * D[i]
  }
  
  # Return derivatives as a concatenated vector
  return(list(c(dNdt, dPdt, dDdt)))
}

#define parameters
u <- 1              # Advection velocity (1/day)
d <- 5              # Diffusion coefficient (m^2/day)
kw <- 0.045         # Light attenuation due to water (m^-1)
kp <- 15e-12         # Light attenuation due to phytoplankton (m^2/cell)
kp_w <- 0           # Light attenuation without phytoplankton (m^2/cell)
Hl <- 20            # Half-saturation constant for light (micromol photons/m^2)
Hn <- 0.1           # Half-saturation constant for nutrients (mmol nutrient/m^3)
gmax <- 1           # Maximum growth rate (1/day)
m <- 0.1            # Mortality rate (1/day)
I0 <- 200           # Surface light intensity (micromol photons/m^2 s)
a <- 1e-9           # Phytoplankton nutrient content (mmol nutrient/cell) 
e <- 0.5            # Nutrient recycling coefficient (#)
Nb <- 350           # Bottom nutrient concentration (mmol nutrient/m^3)
t <- 0.1            # remineralization rate (1/day)

dZ <- 1             # Spatial step size (m)
depth <- 100        # Total depth of the water column (m)

#define grid
Z <- seq(0.5 * dZ, depth - 0.5 * dZ, by = dZ)

#initial conditions
N <- rep(0.5, length(Z))  # Initialize nutrient concentration to zero
P <- rep(0, length(Z))  # Initialize phytoplankton concentration to zero
D <- rep(0, length(Z))  # Initialize detritus concentration to zero
P[10:15] <- 10          # Set initial phytoplankton concentration at specific points
Y <- c(N, P, D)           # Concatenate N and P into a single vector

# Define the time steps
times <- seq(0, 1000, 1)  # Smaller time step for stability

# Solve the ODE
out1 <- ode(y = Y, times = times, func = derivatives, parms = NULL)






              ###Plotting###

# Extract the results from the ODE solution
N_results <- out1[, 2:(length(Z) + 1)]                  # Nutrient concentrations
P_results <- out1[, (length(Z) + 2):(2*length(Z) + 1)]  # Phytoplankton cells
D_results <- out1[, (2*length(Z) + 2):(3*length(Z) + 1)]# Detritus cells

# Create a data frame for plotting
time <- out1[, 1]  # Time steps
depth <- Z         # Depth steps
dataN <- melt(N_results, varnames = c("Time", "Depth"), value.name = "Concentration")
dataP <- melt(P_results, varnames = c("Time", "Depth"), value.name = "Cells")
dataD <- melt(D_results, varnames = c("Time", "Depth"), value.name = "Cells")

# Plot N
ggplot(dataN, aes(x = Time, y = Depth, fill = Concentration)) +
  geom_tile() +  # Create heatmap
  scale_fill_viridis_c(name = "Nutrient\nConcentration [mmol/m^3]") +  # Use viridis color scale
  scale_y_reverse() +  # Reverse y-axis to show depth increasing downward
  labs(x = "Time (days)", y = "Depth (m)", title = "Nutrient Concentration Over Time and Depth") 

# Plot P
ggplot(dataP, aes(x = Time, y = Depth, fill = Cells)) +
  geom_tile() +  # Create heatmap
  scale_fill_viridis_c(name = "Phytoplankton\nConcentration [cells/m^3]") +  # Use viridis color scale
  scale_y_reverse() +  # Reverse y-axis to show depth increasing downward
  labs(x = "Time (days)", y = "Depth (m)", title = "Phytoplankton Concentration Over Time and Depth")

# Plot D
ggplot(dataD, aes(x = Time, y = Depth, fill = Cells)) +
  geom_tile() +  # Create heatmap
  scale_fill_viridis_c(name = "Detritus\nConcentration [cells/m^3]") +  # Use viridis color scale
  scale_y_reverse() +  # Reverse y-axis to show depth increasing downward
  labs(x = "Time (days)", y = "Depth (m)", title = "Detritus Concentration Over Time and Depth")




              ###Plotting last time step###

# Extract the last timestep values
last_time <- tail(time, 1)  # Last timestep
last_N <- N_results[nrow(N_results), ]  # Nutrient concentration at last timestep
last_P <- P_results[nrow(P_results), ]  # Phytoplankton concentration at last timestep
last_D <- D_results[nrow(D_results), ]  # Detritus concentration at last timestep

# Calculate light intensity at the last timestep
last_I <- light(last_P, last_D, Z, I0, kw, kp, dZ)

# Create a data frame for the last timestep
last_timestep_data <- data.frame(
  Depth = -Z,
  Nutrient = last_N,
  Phytoplankton = last_P,
  Detritus = last_D,
  Light = last_I
)

# Reshape the data for plotting
last_timestep_data_long <- melt(last_timestep_data, id.vars = "Depth", variable.name = "Variable", value.name = "Concentration")

# Set up a 1x3 layout for the plots
par(mfrow = c(1, 4))  # 1 row, 4 columns

# Plot 1: Nutrient Concentration
plot(last_timestep_data$Nutrient, last_timestep_data$Depth, 
     type = "l", col = "blue", lwd = 2,
     xlab = "Nutrient Concentration [mmol/m^3]", 
     ylab = "Depth (m)", 
     main = "Nutrient Concentration",
     ylim = c(-100, 0))  # Set y-axis limits

# Plot 2: Phytoplankton Concentration
plot(last_timestep_data$Phytoplankton, last_timestep_data$Depth, 
     type = "l", col = "green", lwd = 2,
     xlab = "Phytoplankton Concentration [cells/m^3]", 
     ylab = "Depth (m)", 
     main = "Phytoplankton Concentration",
     ylim = c(-100, 0))  # Set y-axis limits

# Plot 3: Detritus Intensity
plot(last_timestep_data$Detritus, last_timestep_data$Depth, 
     type = "l", col = "black", lwd = 2,
     xlab = "Detritus Concentration [cells/m^3]", 
     ylab = "Depth (m)", 
     main = "Detritus Concentration",
     ylim = c(-100, 0))  # Set y-axis limits

# Plot 4: Light Intensity
plot(last_timestep_data$Light, last_timestep_data$Depth, 
     type = "l", col = "orange", lwd = 2,
     xlab = "Light Intensity [µmol photons/m^2/s]", 
     ylab = "Depth (m)", 
     main = "Light Intensity",
     ylim = c(-100, 0))  # Set y-axis limits
# Reset the plotting layout to default
par(mfrow = c(1, 1))