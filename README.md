# Computational-marine-ecological-modelling
library(deSolve)
library(ggplot2)
library(reshape2)

#define parameters
u <- 0.04           # Advection velocity
D <- 5              # Diffusion coefficient
kw <- 0.2           # Light attenuation due to water
kp <- 15 * 10^(-12) # Light attenuation due to phytoplankton
H <- 30             # Half-saturation constant for light
gmax <- 0.04*24       # Maximum growth rate
m <- 0.01*24           # Mortality rate
I0 <- 350           # Surface light intensity

dZ <- 1# Spatial step size
depth <- 120  # Total depth of the water column

#define grid
Z <- seq(0.5 * dZ, depth - 0.5 * dZ, by = dZ)

#initial conditions
P <- rep(0, length(Z))  # Initialize concentration to zero
P[50:59] <- 1  # Set initial concentration at specific points
P0 <- P

# Define the light function
light <- function(P, Z, I0, kw, kp, dZ) {
  L <- cumsum(kp * P * dZ)  # Cumulative light attenuation due to phytoplankton
  I <- I0 * exp(-kw*Z - L)  # Light intensity at depth Z
  return(I)
}

# Define the derivatives function
derivatives <- function(t, P, parms) {
  # Calculate light intensity
  I <- light(P, Z, I0, kw, kp, dZ)
  
  # Define boundaries
  Ja <- rep(0, length(Z) + 1)
  Jd <- rep(0, length(Z) + 1)
  
  # Define fluxes
  for (i in 2:length(Z)) {
    Ja[i] <- u * P[i - 1]
    Jd[i] <- -D * (P[i] - P[i - 1]) / dZ
  }
  
  # Boundary conditions (no flux at top and bottom)
  Ja[1] <- 0
  Jd[1] <- 0
  Ja[length(Z) + 1] <- 0
  Jd[length(Z) + 1] <- 0
  
  J <- Ja + Jd
  
  dPdt <- rep(0, length(Z))
  
  for (i in 1:length(Z)) {
    dPdt[i] <- (-(J[i + 1] - J[i]) / dZ) + (-m + gmax * I[i] / (H + I[i])) * P[i]
  }
  
  return(list(dPdt))
}

# Define the time steps
times <- seq(0, 200, 0.5)  # Smaller time step for stability

# Solve the ODE
out1 <- ode(y = P, times = times, func = derivatives, parms = NULL)





##Plotting
# Convert the output to a data frame
out1_df <- as.data.frame(out1)
colnames(out1_df) <- c("time", paste("Z", 1:length(Z), sep = "_"))

# Reshape the data for plotting
out1_melt <- melt(out1_df, id.vars = "time", variable.name = "depth", value.name = "P")

# Convert depth to numeric
out1_melt$depth <- as.numeric(gsub("Z_", "", out1_melt$depth)) * dZ

# Define custom color scale
custom_colors <- c("darkblue", "yellow", "green3")  # Gradient from blue to yellow to green

# Plot the results with custom color scale
ggplot(out1_melt, aes(x = time, y = depth, fill = P)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradientn(
    colors = custom_colors,  # Gradient colors
    values = scales::rescale(c(0, 0.5, 1)),  # Rescale to ensure smooth transition
    guide = "colourbar"  # Add a color legend
  ) +
  labs(x = "Time", y = "Depth", fill = "Concentration (P)") +
  theme_minimal()
