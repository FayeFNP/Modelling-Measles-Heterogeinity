### The Setup --------------------------------------------------------------------

# install and load necessary packages
if(!require(viridis)){install.packages('viridis')}
library(viridis)
if(!require(devtools)){install.packages('devtools')}
library(devtools)
if(!require(tidyverse)){install.packages('tidyverse')}
library(tidyverse)
if(!require(parallel)){install.packages('parallel')}
library(parallel)
if(!require(tidyr)){install.packages('tidyr')}
library(tidyr)
if(!require(gganimate)){install.packages('gganimate')}
library(gganimate)
if(!require(patchwork)){install.packages('patchwork')}
library(patchwork)
if(!require(gifski)){install.packages('gifski')}
library(gifski)
if(!require(scales)){install.packages('scales')}
library(scales)

# Load custom model functions and data
#load_all("./SEITmeta")
require(SEITmeta)
load("../Data/MeaslesPosterior.RData")
source('../Rinclude/Model_Functions.R')

# Model parameters
TE = 8
kE = 8
TI = 5
kI = 5
TT = Inf

# Distance matrix to calculate geographic extent
dm = as.matrix(dist(z$xy))

# For initial conditions, setup up elimination state coverage of 97% (higher rate of external imports)
# Start with no exposed or infectious individuals, sparking of outbreaks is only influenced by imported cases

init <- vacc_scenario_init(0.97)
init$T = init$T + init$E + init$I
init$E = 0
init$I = 0
save(init, file = "init.RData")
# load("init.RData")


### Swansea Region Homogeneous Scenarios ------------------------------------------
n_reps <- 10
v_grid <- seq(0.85, 0.99, by = 0.02)
import_grid <- 10^seq(-8, -3, by = 1)
pars <- expand.grid(v = v_grid, import_rate = import_grid)

# Swansea region definition == 29 closest districts surrounding Swansea (Swansea = district 959)
SwanseaR <- sort.int(dm[959, ], index.return = TRUE)$ix[1:29]

# To control for multiple seeds sparking outbreaks everywhere, seeding in surrounding Swansea districts
seed_swanseaR <- function(import_rate) {
  import_vec = numeric(length(init$S))
  import_vec[SwanseaR] = import_rate  
  return(import_vec)
}

# Run simulations with 10 replicates per scenario
swansea_outputR <- mclapply(1:nrow(pars), function(i) {
  replicate(n_reps, {
    vacc_scenario_ext(rep(pars$v[i],982), init, import_rate = seed_swanseaR(pars$import_rate[i]))
  }, simplify = FALSE) %>% bind_rows()
}, mc.cores = 32)

# Combine results 
swansea_outputR <- bind_rows(swansea_outputR, .id = "param_combo") %>%
  mutate(param_combo = as.integer(param_combo)) %>%
  left_join(pars %>% mutate(param_combo = row_number()), by = "param_combo") %>%
  relocate(v, import_rate, .before = everything())

save(swansea_outputR, file = "swansearegion_homogeneous_results.RData")
load("swansearegion_homogeneous_results.RData")


# Plot Affected Districts
ggplot(swansea_outputR, aes(x = v, y = prevalence_infected_districts, colour = factor(import_rate), fill = factor(import_rate))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.4, colour = NA) +
  facet_wrap(~ import_rate, ncol = 3) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Prevalence of Affected Districts by Vaccine Coverage",
    x = "Vaccine Coverage",
    y = "No. of Affected Districts",
    colour = "Import Rate",
    fill = "Import Rate"
  ) +
  theme_minimal()


# plot max spatial extent
swansea_outputR_clean <- swansea_outputR %>%
  filter(is.finite(maximum_extent))

ggplot(swansea_outputR_clean, aes(x = v, y = maximum_extent, colour = factor(import_rate), fill = factor(import_rate))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.4, colour = NA) +
  facet_wrap(~ import_rate, ncol = 3) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Geographic Spatial Extent of Outbreaks by Vaccine Coverage",
    x = "Vaccine Coverage",
    y = "Max Distance Between Affected Districts (km)",
    colour = "Import Rate",
    fill = "Import Rate"
  ) +
  theme_minimal()

# Plot peak size 
swansea_outputR_clean2 <- swansea_outputR %>% 
  filter(is.finite(time_to_first_outbreak), is.finite(peak_first_outbreak))

ggplot(swansea_outputR_clean2, aes(x = v, y = peak_first_outbreak, colour = factor(import_rate), fill = factor(import_rate))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.4, colour = NA) +
  facet_wrap(~ import_rate, ncol = 3) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Size of First Outbreak across Vaccine Coverage",
    x = "Vaccine Coverage",
    y = "No of Cases at Peak",
    colour = "Import Rate",
    fill = "Import Rate"
  ) +
  theme_minimal()


# Plot time taken
ggplot(swansea_outputR_clean2, aes(x = v, y = time_to_first_outbreak, colour = factor(import_rate), fill = factor(import_rate))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.4, colour = NA) +
  facet_wrap(~ import_rate, ncol = 3) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Time to First Outbreak across Vaccine Coverage",
    x = "Vaccine Coverage",
    y = "No of Weeks",
    colour = "Import Rate",
    fill = "Import Rate"
  ) +
  theme_minimal()



### Swansea Region Heterogeneity Simulations -----------------------------------------------
set.seed(432)
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 100                        

hetero_swR <- mclapply(alphas, function(a) {
  replicate(n_reps, {
    pattern <- rep(TRUE, 982)
    pattern[sort.int(dm[959,],index=T)$ix[1:29]] = FALSE
    vacc_coverage <- mk_vacc_coverage(0.97, pattern, a)
    
    # Store vh and vl for analysis
    s <- sum(pattern)/length(pattern)
    vh <- (0.97 * (1 + 2 * a * (1 - s)))    # for pattern == TRUE
    vl <- (0.97 * (1 - 2 * a * s))          # for pattern == FALSE
    
    # To clarify, low coverage is in Swansea region (FALSE) + seeding also in the region 
    import_rate <- rep(0,982)
    import_rate[sort.int(dm[959,],index=T)$ix[1:29]] = (10^-6)
    res <- vacc_scenario_ext(vacc_coverage, init, import_rate)
    res$alpha <- a
    res$vh <- vh
    res$vl <- vl
    return(res)
  }, simplify = FALSE) %>% bind_rows()
}, mc.cores = 32) %>% bind_rows()

save(hetero_swR, file = "swanseaR_results.RData")
load("swanseaR_results.RData")

# Summarise with SD
hetero_summary_swR <- hetero_swR %>%
  group_by(alpha) %>%
  summarise(
     mean_prevalence = mean(prevalence_infected_districts),
    sd_prevalence = sd(prevalence_infected_districts),
    
    mean_extent = mean(maximum_extent),
    sd_extent = sd(maximum_extent),
    
    mean_peak_size = mean(peak_first_outbreak),
    sd_peak_size = sd(peak_first_outbreak),
    
    mean_time = mean(time_to_first_outbreak),
    sd_time = sd(time_to_first_outbreak),
    
    .groups = "drop"
  )

save(hetero_summary_swR, file = "swRhetero_summary_sd.RData")
load("swRhetero_summary_sd.RData")


# Prevalence of affected districts
p4 <- ggplot(hetero_summary_swR, aes(x = alpha, y = mean_prevalence)) +
  geom_ribbon(aes(ymin = pmax(mean_prevalence - sd_prevalence, 0),
                  ymax = mean_prevalence + sd_prevalence),
              fill = "grey80", alpha = 0.4) +
  geom_line() +
  geom_point() +
  labs(
    title = "D. Prevalence of Affected Districts",
    x = expression("Vaccine Heterogeneity"~(alpha)),
    y = "No. of Affected Districts"
  ) +
  theme_minimal()

# Spatial extent
p3 <- ggplot(hetero_summary_swR, aes(x = alpha, y = mean_extent)) +
  geom_ribbon(aes(ymin = pmax(mean_extent - sd_extent, 0),
                  ymax = mean_extent + sd_extent),
              fill = "grey80", alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    title = "C. Spatial Extent of Outbreaks",
    x = expression("Vaccine Heterogeneity"~(alpha)),
    y = "Max Distance Between Infected Districts (km)"
  ) +
  theme_minimal()

# First outbreak size
p2 <- ggplot(hetero_summary_swR, aes(x = alpha, y = mean_peak_size)) +
  geom_ribbon(aes(ymin = pmax(mean_peak_size - sd_peak_size, 0),
                  ymax = mean_peak_size + sd_peak_size),
              fill = "grey80", alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    title = "B. Size of First Outbreak",
    x = expression("Vaccine Heterogeneity"~(alpha)),
    y = "No. of Cases at Peak"
  ) +
  theme_minimal()

# Time to first outbreak
p1 <- ggplot(hetero_summary_swR, aes(x = alpha, y = mean_time)) +
  geom_ribbon(aes(ymin = pmax(mean_time - sd_time, 0),
                  ymax = mean_time + sd_time),
              fill = "grey80", alpha = 0.4) +
  geom_point() +
  geom_line() +
  labs(
    title = "A. Time to First Outbreak",
    x = expression("Vaccine Heterogeneity"~(alpha)),
    y = "No. of Weeks"
  ) +
  theme_minimal()

# Combine into a 2x2 layout
combined_plot2 <- (p1 | p2) / (p3 | p4)

# Display the figure
print(combined_plot2)

# Save to file
ggsave("swheterogeneity_4panel_plot.pdf", combined_plot2, width = 12, height = 10)


### Swansea Region Maps --------------------------------------------------------
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 100

# summarise cases per district
summarise_cases <- function(alpha) {
    pattern <- rep(TRUE, 982)
    pattern[sort.int(dm[959,],index=T)$ix[1:29]] = FALSE
    vacc_coverage <- mk_vacc_coverage(0.97, pattern, alpha)
    import_rate <- rep(0,982)
    import_rate[sort.int(dm[959,],index=T)$ix[1:29]] = (10^-6)
  
replicate(n_reps, {
    # Run the model
    cases <- vacc_scenario_extC(vacc_coverage, init, import_rate)
    # Sum across rows (i.e., over time) to get cumulative cases per district
    cumulative_cases <- colSums(cases)
    # Create a tibble for mapping
    tibble(
      district = 1:length(cumulative_cases),
      cumulative_cases = cumulative_cases,
      alpha = alpha
    )
  }, simplify = FALSE) %>% bind_rows()
}

# Run in parallel
map_cases_sw <- mclapply(alphas, summarise_cases, mc.cores = 6) %>%
  bind_rows()

save(map_cases_sw, file = "Swansea_maps.RData")

# Merge with spatial coordinates
Swansea_map_data <- map_cases_sw %>% 
  left_join(
    as.data.frame(z$xy) %>% 
      setNames(c("x", "y")) %>% 
      mutate(district = 1:nrow(z$xy)),
    by = "district"
  )

save(Swansea_map_data, file = "Swansea_map_data.RData")

load("Swansea_map_data.RData")
# Plot faceted maps
ggplot(Swansea_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
  geom_point(size = 1) +
  scale_colour_viridis_c(option = "plasma", trans = "log1p", breaks = c(0, 5000,15000),
                         labels = c("0", "5000","15000")) +
  facet_wrap(~ alpha) +
  coord_equal() +
  labs(
    title = "Cumulative Cases Across Districts by Alpha",
    colour = "Cases"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


#Plot an animated map of cumulative cases over time
Swansea_point_map <- ggplot(Swansea_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
  geom_point(size = 1) +
  scale_colour_viridis_c(option = "plasma", trans = "log1p", breaks = c(0, 5000,10000,15000,20000),
                         labels = c("0", "5000", "10000", "15000", "20000")) +
  coord_equal() +
  labs(
    title = "Cumulative Cases Across Districts (Alpha = {closest_state})",
    colour = "Cases"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  transition_states(alpha, transition_length = 2, state_length = 1) +
  ease_aes("cubic-in-out")


casemap_sw <- animate(Swansea_point_map, nframes = 100, fps = 10, width = 800, height = 600, 
renderer = gifski_renderer())
anim_save("Swansea_cumulative_cases.gif", casemap_sw)


### Swansea Region Time Series ------------------------------------------------------------------
# Temporal trajectories of key scenarios
# Define alpha values
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 10

# Spatial Ring Aggregates

# Get coordinates of Swansea centroid (district 959)
swansea_coords <- z$xy[959, ]
coords_df <- as.data.frame(z$xy)
names(coords_df) <- c("x", "y")
coords_df$district <- 1:nrow(coords_df)

# Compute distance from Swansea centroid
coords_df <- coords_df %>%
  mutate(distance = sqrt((x - swansea_coords[1])^2 + (y - swansea_coords[2])^2),
         zone = case_when(
           district == 959 ~ "Swansea",
           distance <= 100 ~ "Nearby (≤100km)",
           TRUE ~ "Distant (>100km)"
         ))

ts_all_alpha_zones <- mclapply(alphas, function(a) {
  bind_rows(lapply(seq_len(n_reps), function(r) {
    pattern <- rep(TRUE,982)
    pattern[sort.int(dm[959,],index=T)$ix[1:29]] = FALSE
    vacc_coverage <- mk_vacc_coverage(0.97, pattern, a)
    import_rate <- rep(0,982)
    import_rate[sort.int(dm[959,],index=T)$ix[1:29]] = (10^-6)
    
    ts_output <- vacc_scenario_extT(vacc_coverage, init, import_rate)
    
    n_patches <- 982
    n_weeks <- length(unique(ts_output$time))
    if (!"district" %in% names(ts_output)) {
      ts_output$district <- rep(1:n_patches, times = n_weeks)
    }
    
    ts_output$zone <- coords_df$zone[ts_output$district]
    ts_output$replicate <- r
    
    ts_output %>%
      group_by(time, zone, replicate) %>%
      summarise(total_cases = sum(cases), .groups = "drop") %>%
      mutate(alpha = a)
  }))
}, mc.cores = 6) %>% bind_rows()


ggplot(ts_all_alpha_zones, aes(x = time, y = total_cases, colour = zone, group = interaction(zone, replicate))) +
  geom_path(alpha = 0.1, linewidth = 0.4) +
  facet_wrap(~ alpha, ncol = 4) +
  labs(
    title = "Epidemic Trajectories by Zone Across Increasing Spatial Heterogeneity",
    subtitle = "Mean cumulative cases per zone",
    x = "Time (years)",
    y = "Cumulative Cases",
    colour = "Zone"
  ) +
  scale_x_continuous(
    breaks = seq(104, 520, by = 104),  # Tick every 2 years (52 weeks × 2)
    labels = function(x) paste0(x / 52)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 3)
  ) +
  scale_colour_manual(
    values = c(
      "Swansea" = "red",
      "Nearby (≤100km)" = "orange",
      "Distant (>100km)" = "forestgreen"  
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )