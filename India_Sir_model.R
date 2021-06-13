### SIR model. THis is using India incidences for Covid-19. Input required is 
# cumulative incidence and total population. It calculates beta and gamma using ODE.
# it also calculate reproduction rate R0 that means how many people get affected from one person.
# the code also gives peak date and maximum people that will get infected. 
#setwd("C:\\Users\\Talal Mohd\\Desktop\\tayyab\\TERI\\Internship(Minor Project)\\SIR\\India")
India <- read.table('India_nation_level_daily.csv', header=TRUE, sep=',')
India$Date <- as.Date(India$Date,format = "%m/%d/%y")
N <- 1387297452 # total population
SIR <- function(time, country, parameters) {
  par <- as.list(c(country, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}
library(deSolve)
library(lubridate)
sir_start_date <- "2020-01-30"
sir_end_date <- "2020-06-29"                              

Infected <- India$I
Day <- 1:(length(Infected))

# now specify initial values for S, I and R

init <- c(
  S = N - Infected[1],
  I = Infected[1],
  R = 0
)

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[, 3]
  sum((Infected - fit)^2)
}

# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0


Opt <- optim(c(0.1, 0.1),
             RSS,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)


# check for convergence
Opt$message
#Convergence is confirmed. Now we can examine the fitted values for ββ and γγ:

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date))
fitted_cumulative_incidence <- data.frame(ode(
  y = init, times = t,
  func = SIR, parms = Opt_par
))

library(dplyr)
fitted_cumulative_incidence <- fitted_cumulative_incidence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t - 1),
    cumulative_incident_cases = Infected
  )

# plot the data
library(ggplot2)
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = cumulative_incident_cases), colour = "blue") +
  labs(
    y = "Cumulative incidence",
    title = "COVID-19 fitted vs observed cumulative incidence, India",
    subtitle = "(Red = fitted from SIR model, blue = observed)"
  ) +
  theme_minimal()

R0 <- as.numeric(Opt_par[1] / Opt_par[2])
# plot in log scale
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = cumulative_incident_cases), colour = "blue") +
  labs(
    y = "Cumulative incidence",
    title = "COVID-19 fitted vs observed cumulative incidence, India(Log Scale)",
    subtitle = "(Red = fitted from SIR model, blue = observed)"
  ) +
  theme_minimal() +
  scale_y_log10(labels = scales::comma)


# time in days for predictions
t <- 1:270

# get the fitted values from our SIR model
fitted_cumulative_incidence <- data.frame(ode(
  y = init, times = t,
  func = SIR, parms = Opt_par
))


# add a Date column and join the observed incidence data
fitted_cumulative_incidence <- fitted_cumulative_incidence %>%
  mutate(
    Date = ymd(sir_start_date) + days(t - 1),
    country = "India",
    cumulative_incident_cases = c(Infected, rep(NA, length(t) - length(Infected)))
  )

# plot the data
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I, colour = "red")) +
  geom_line(aes(y = S, colour = "black")) +
  geom_line(aes(y = R, colour = "green")) +
  geom_point(aes(y = cumulative_incident_cases, colour = "blue")) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = "Persons", title = "COVID-19 fitted vs observed cumulative incidence, India") +
  scale_colour_manual(name = "", values = c(red = "red", black = "black", green = "green", blue = "blue"
  ), labels = c("Susceptible","Recovered", "Observed", "Infectious" )) +
  theme_minimal()


#The same graph in log scale for the y-axis and with a legend for better readability:
# plot the data
fitted_cumulative_incidence %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I, colour = "red")) +
  geom_line(aes(y = S, colour = "black")) +
  geom_line(aes(y = R, colour = "green")) +
  geom_point(aes(y = cumulative_incident_cases, colour = "blue")) +
  scale_y_log10(labels = scales::comma) +
  labs(
    y = "Persons",
    title = "COVID-19 fitted vs observed cumulative incidence, India(Log Scale)"
  ) +
  scale_colour_manual(
    name = "",
    values = c(red = "red", black = "black", green = "green", blue = "blue"),
    labels = c("Susceptible", "Observed", "Recovered", "Infectious")
  ) +
  theme_minimal()

#Other interesting statistics can be computed from the fit of our model. For example:
# the peak of the pandemic
#	the number of severe cases
#	the number of people in need of intensive care
#	the number of deaths
fit <- fitted_cumulative_incidence

# peak of pandemic
fit[fit$I == max(fit$I), c("Date", "I")]

# severe cases
max_infected <- max(fit$I)
# cases with need for intensive care
intensive_care_patient <- max_infected * 0.02 # assuming 2% infected people will need intensive care
# if we consider death rate as total death/total infected till June 29
death_rate = India$Deaths[dim(India)[1]]/India$Confirmed[dim(India)[1]]
total_death = max_infected * death_rate
