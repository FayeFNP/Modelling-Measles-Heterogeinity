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
#load("init.RData")

### London Homogeneous Scenarios --------------------------------------------------------
n_reps <- 10
v_grid <- seq(0.85, 0.99, by = 0.02)
import_grid <- 10^seq(-8, -3, by = 1)
pars <- expand.grid(v = v_grid, import_rate = import_grid)

# To control for multiple seeds sparking outbreaks everywhere, seeding only in London
seed_london <- function(import_rate) {
  import_vec = numeric(length(init$S))
  import_vec[1:29] = import_rate  
  return(import_vec)
}

# Run simulations with 10 replicates per scenario
london_output <- mclapply(1:nrow(pars), function(i) {
  replicate(n_reps, {
    vacc_scenario_ext(rep(pars$v[i],982), init, import_rate = seed_london(pars$import_rate[i]))
  }, simplify = FALSE) %>% bind_rows()
}, mc.cores = 32)

# Combine results 
london_output <- bind_rows(london_output, .id = "param_combo") %>%
  mutate(param_combo = as.integer(param_combo)) %>%
  left_join(pars %>% mutate(param_combo = row_number()), by = "param_combo") %>%
  relocate(v, import_rate, .before = everything())

save(london_output, file = "london_homogeneous_results.RData")
load("london_homogeneous_results.RData")


# Plot Affected Districts
ggplot(london_output, aes(x = v, y = prevalence_infected_districts, colour = factor(import_rate), fill = factor(import_rate))) +
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
london_output_clean <- london_output %>%
  filter(is.finite(maximum_extent))

ggplot(london_output_clean, aes(x = v, y = maximum_extent, colour = factor(import_rate), fill = factor(import_rate))) +
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
london_output_clean2 <- london_output %>% 
  filter(is.finite(time_to_first_outbreak), is.finite(peak_first_outbreak))

ggplot(london_output, aes(x = v, y = peak_first_outbreak, colour = factor(import_rate), fill = factor(import_rate))) +
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
ggplot(london_output_clean2, aes(x = v, y = time_to_first_outbreak, colour = factor(import_rate), fill = factor(import_rate))) +
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


### London Heterogeneity Simulations ----------------------------------------------------

set.seed(432)
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 100                        

hetero_ld <- mclapply(alphas, function(a) {
  replicate(n_reps, {
    pattern <- rep(TRUE, 982)
    pattern[1:29] <- FALSE
    vacc_coverage <- mk_vacc_coverage(0.97, pattern, a)
    
    # Store vh and vl for analysis
    s <- sum(pattern)/length(pattern)
    vh <- (0.97 * (1 + 2 * a * (1 - s)))    # for pattern == TRUE
    vl <- (0.97 * (1 - 2 * a * s))          # for pattern == FALSE
    
    # To clarify, low coverage is in London (FALSE) + seeding also in London 
    import_rate <- seed_london(10^-6)
    res <- vacc_scenario_ext(vacc_coverage, init, import_rate)
    res$alpha <- a
    res$vh <- vh
    res$vl <- vl
    return(res)
  }, simplify = FALSE) %>% bind_rows()
}, mc.cores = 32) %>% bind_rows()

save(hetero_ld, file = "london_results.RData")
load("london_results.RData")

# Summarise with SD
hetero_summary_ld <- hetero_ld %>%
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

save(hetero_summary_ld, file = "hetero_summary_sd.RData")
load("hetero_summary_sd.RData")

# Prevalence of affected districts
p4 <- ggplot(hetero_summary_ld, aes(x = alpha, y = mean_prevalence)) +
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
p3 <- ggplot(hetero_summary_ld, aes(x = alpha, y = mean_extent)) +
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
p2 <- ggplot(hetero_summary_ld, aes(x = alpha, y = mean_peak_size)) +
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
p1 <- ggplot(hetero_summary_ld, aes(x = alpha, y = mean_time)) +
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

# Combine all plots into a 2x2 layout
combined_plot <- (p1 | p2) / (p3 | p4)

# Display the combined figure
print(combined_plot)

# Save to file
ggsave("heterogeneity_4panel_plot.pdf", combined_plot, width = 12, height = 10)


### London Maps ---------------------------------------------------------------------
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 100

# summarise cases per district
summarise_cases <- function(alpha) {
  pattern <- c(rep(FALSE, 29), rep(TRUE, 982 - 29))
  vacc_coverage <- mk_vacc_coverage(0.97, pattern, alpha)
  import_rate <- seed_london(10^-6)
  
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
map_cases_ld <- mclapply(alphas, summarise_cases, mc.cores = 6) %>%
  bind_rows()

save(map_cases_ld, file = "London_maps.RData")

# Merge with spatial coordinates
ld_map_data <- map_cases_ld %>% 
  left_join(
    as.data.frame(z$xy) %>% 
      setNames(c("x", "y")) %>% 
      mutate(district = 1:nrow(z$xy)),
    by = "district"
  )

save(ld_map_data, file = "ld_map_data.RData")
load("ld_map_data.RData")
# Plot faceted maps
ggplot(ld_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
  geom_point(size = 1) +
  scale_colour_viridis_c(option = "plasma", trans = "log1p", breaks = c(0, 5000,20000),
                         labels = c("0", "5000", "20000")) +
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
ld_point_map <- ggplot(ld_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
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

casemap_ld <- animate(ld_point_map, nframes = 100, fps = 10, width = 800, height = 600, 
renderer = gifski_renderer())
anim_save("London_cumulative_cases.gif", casemap_ld)


### Location Time Series -------------------------------------------------------
alphas <- seq(0, 0.3, by = 0.05)
n_reps <- 10

# Get coordinates of London centroid (average of districts 1:29)
london_coords <- colMeans(z$xy[1:29, ])
coords_df <- as.data.frame(z$xy)
names(coords_df) <- c("x", "y")
coords_df$district <- 1:nrow(coords_df)

# Compute distance from London centroid
coords_df <- coords_df %>%
  mutate(distance = sqrt((x - london_coords[1])^2 + (y - london_coords[2])^2),
         zone = case_when(
           district %in% 1:29 ~ "London",
           distance <= 100 ~ "Nearby (≤100km)",
           TRUE ~ "Distant (>100km)"
         ))
 
ts_all_alpha_zones <- mclapply(alphas, function(a) {
  bind_rows(lapply(seq_len(n_reps), function(r) {
    pattern <- c(rep(FALSE, 29), rep(TRUE, 982 - 29))
    vacc_coverage <- mk_vacc_coverage(0.97, pattern, a)
    import_rate <- seed_london(1e-6)
    
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
      "London" = "red",
      "Nearby (≤100km)" = "orange",
      "Distant (>100km)" = "forestgreen"  
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )
