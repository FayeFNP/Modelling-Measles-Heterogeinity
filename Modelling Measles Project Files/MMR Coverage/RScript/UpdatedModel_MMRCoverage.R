# The Setup --------------------------------------------------------------------
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
write.csv(z$xy, "XYlatlong.csv", row.names = FALSE)
# Load MMR coverage data with UD-to-UTLA dictionary

load('../Data/DictionaryUD_to_UTLA.RData')
MMRCoverage <- MMR_coverage
#MMRCoverage <- read_csv("../Data/DictionaryUD_to_UTLA.csv")

# Arrange district ID from 1 to 982
MMRCoverage <- MMRCoverage %>%
  arrange(UD) %>%
  mutate(coverage = coverage / 100)  # Convert to proportion

# Coverage vector for 982 districts
MMRcoverage <- numeric(982)
MMRcoverage[MMRCoverage$UD] <- MMRCoverage$coverage

# Model parameters
TE = 8
kE = 8
TI = 5
kI = 5
TT = Inf

# Distance matrix to calculate geographic extent
dm = as.matrix(dist(z$xy))   # 982 x 982 matrix

# For initial conditions, setup up elimination state coverage of 97%
# Start with no exposed or infectious individuals, sparking of outbreaks is only influenced by the imported cases

init_MMR <- vacc_scenario_init(0.97)
init_MMR$T = init_MMR$T + init_MMR$E + init_MMR$I
init_MMR$E = 0
init_MMR$I = 0
save(init_MMR, file = "init_MMR.RData")
#load("init_MMR.RData")

# MMR Simulation ----------------------------------------------------
# Group districts with coverage < 80% outside London 
# Lowest predicted coverage results from the binomial regression model
hotspot_districts <- c(892, 942, 972)

# Seed hotspots
seed_hotspots <- function(import_rate, hotspot_districts) {
  import_vec = numeric(982)
  import_vec[hotspot_districts] <- import_rate  
  return(import_vec)
}

import_rates <- 10^seq(-8, -3, by = 1)
n_reps <- 100
set.seed(432)
                     

hetero_mmr <- mclapply(import_rates, function(imp_rate) {
  replicate(n_reps, {
    vacc_coverage <- MMRcoverage 
    import_vec <- seed_hotspots(imp_rate, hotspot_districts)
    
    res <- vacc_scenario_ext(vacc_coverage, init_MMR, import_vec)
    res$import_rate <- imp_rate
    
    return(res)
  }, simplify = FALSE) %>% bind_rows()
}, mc.cores = 32) %>% bind_rows()

save(hetero_mmr, file = "MMR_results.RData")
#load("MMR_results.RData")

# Summarise with SD
hetero_mmr2 <- hetero_mmr %>%
  group_by(import_rate) %>%
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

save(hetero_mmr2, file = "hetero_sd.RData")
#load("hetero_sd.RData")


# Combine plots into a 2x2 layout using patchwork

# Time to first outbreak
p1 <- ggplot(hetero_mmr2, aes(x = import_rate, y = mean_time)) +
  geom_ribbon(aes(ymin = pmax(0, mean_time - sd_time),
                  ymax = mean_time + sd_time),
              fill = "grey80", alpha = 0.4) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  labs(
    title = "A. Time to First Outbreak",
    x = "Log Import Rate",
    y = "No. of Weeks"
  ) +
  theme_minimal()

# First outbreak size
p2 <- ggplot(hetero_mmr2, aes(x = import_rate, y = mean_peak_size)) +
  geom_ribbon(aes(ymin = pmax(0, mean_peak_size - sd_peak_size),
                  ymax = mean_peak_size + sd_peak_size),
              fill = "grey80", alpha = 0.4) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  labs(
    title = "B. Size of First Outbreak",
    x = "Log Import Rate",
    y = "No. of Cases at Peak"
  ) +
  theme_minimal()

# Prevalence of affected districts
p3 <- ggplot(hetero_mmr2, aes(x = import_rate, y = mean_prevalence)) +
  geom_ribbon(aes(ymin = pmax(0, mean_prevalence - sd_prevalence),
                  ymax = mean_prevalence + sd_prevalence),
              fill = "grey80", alpha = 0.4) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  labs(
    title = "C. Prevalence of Affected Districts",
    x = "Log Import Rate",
    y = "No. of Affected Districts"
  ) +
  theme_minimal()

# Combine into a 3x1 layout (horizontal row)
MMRplot <- p1 | p2 | p3

# Display the combined figure
print(MMRplot)

# Save to file
ggsave("MMR_import_3panel_plot.pdf", MMRplot, width = 14, height = 5)
ggsave("MMR_import_3panel_plot.png", MMRplot, width = 14, height = 5, dpi = 300)

# Time Series ------------------------------------------------------------------
n_reps <- 10
import_rates <- 10^seq(-8, -3, by = 1)

set.seed (567)
ts_imports <- mclapply(import_rates, function(ir) {
  replicate(n_reps, {
    vacc_coverage <- MMRcoverage
    import_rate <- seed_hotspots(ir, hotspot_districts)  
    ts_output <- vacc_scenario_extT(vacc_coverage, init_MMR, import_rate)
    ts_summary <- ts_output %>%
      group_by(time) %>%
      summarise(total_cases = sum(cases), .groups = "drop") %>%
      mutate(import_rate = ir)
    return(ts_summary)
  }, simplify = FALSE) %>% bind_rows()
}) %>% bind_rows()

#Plot animated time series for each alpha
p <- ggplot(ts_imports, aes(x = time, y = total_cases, group = import_rate)) +
  geom_line(color = "steelblue") +
  labs(
    title = "Epidemic Trajectory Across Varying Import Rates",
    subtitle = "Import Rate = {closest_state}",
    x = "Time (weeks)",
    y = "Total Cases"
  ) +
  theme_minimal() +
  transition_states(import_rate, transition_length = 2, state_length = 1) +
  ease_aes("cubic-in-out")

anim <- animate(p, nframes = 100, fps = 10, width = 800, height = 600, 
renderer = gifski_renderer())
anim_save("epidemic_trajectory.gif", anim)

# Faceted time series plot 
ts_imports$import_rate_label <- factor(
  paste0("Import rate = ", formatC(ts_imports$import_rate, format = "e", digits = 1)),
  levels = paste0("Import rate = ", formatC(sort(unique(ts_imports$import_rate)), format = "e", digits = 1))
)

ggplot(ts_imports, aes(x = time, y = total_cases, group = interaction(import_rate_label), colour = import_rate_label)) +
  geom_line(alpha = 0.3, linewidth = 0.4, show.legend = FALSE) +
  facet_wrap(~ import_rate_label, ncol = 3) +
  labs(
    title = "Epidemic Trajectories under Real MMR Coverage",
    subtitle = "Varying Importation Rates",
    x = "Time (weeks)",
    y = "Total Cases"
  ) +
  theme_minimal()

ggsave("facet_import_trajectories.png", width = 12, height = 8, dpi = 300)


# MMR Maps ---------------------------------------------------------------------

import_rate <- 10^seq(-8, -3, by = 1)
n_reps <- 100
hotspot_districts <- c(892, 942, 972) 

# seed only hotspots
seed_hotspots <- function(import_rate, hotspot_districts) {
  import_vec = numeric(982)
  import_vec[hotspot_districts] <- import_rate 
  return(import_vec)
}

# summarise cases per district
summarise_MMRcases <- function(import_rate) {
  import_vec <- seed_hotspots(import_rate, hotspot_districts)
  replicate(n_reps, {
    # Run the model
    cases <- vacc_scenario_extC(MMRcoverage, init_MMR, import_vec)
    # Sum over time
    cumulative_cases <- colSums(cases)
    # Create a tibble for mapping
    tibble(
      district = 1:length(cumulative_cases),
      cumulative_cases = cumulative_cases,
      import_rate = import_rate
    )
  }, simplify = FALSE) %>% bind_rows()
}

set.seed (123)
# Run in parallel
map_MMRcases <- mclapply(import_rates, summarise_MMRcases, mc.cores = 6) %>%
  bind_rows()

save(map_MMRcases, file = "MMR_maps.RData")
#load("MMR_maps.RData")

# Merge with spatial coordinates
MMR_map_data <- map_MMRcases %>% 
  left_join(
    as.data.frame(z$xy) %>% 
      setNames(c("x", "y")) %>% 
      mutate(district = 1:nrow(z$xy)),
    by = "district"
  )

# Add log-scaled import label
MMR_map_data$import_rate_label <- factor(
  paste0("Import = ", formatC(MMR_map_data$import_rate, format = "e", digits = 1)),
  levels = paste0("Import = ", formatC(sort(unique(MMR_map_data$import_rate)), format = "e", digits = 1))
)

save(MMR_map_data, file = "MMR_map_data.RData")

# Plot faceted maps
ggplot(MMR_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
  geom_point(size = 1) +
  scale_colour_viridis_c(option = "plasma", trans = "log1p", breaks = c(0, 10000,20000,50000,100000),
                         labels = c("0", "10000", "20000", "50000", "100000")) +
  facet_wrap(~ import_rate_label) +
  coord_equal() +
  labs(
    title = "Cumulative Cases Across Districts",
    colour = "Cases"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


#Plot an animated map of cumulative cases over time
MMR_point_map <- ggplot(MMR_map_data, aes(x = x, y = y, colour = cumulative_cases)) +
  geom_point(size = 1) +
  scale_colour_viridis_c(option = "plasma", trans = "log1p", breaks = c(0, 10000,100000,1000000),
                         labels = c("0", "10K", "100K", "1M")) +
  coord_equal() +
  labs(
    title = "Cumulative Cases Across Districts (Import = {closest_state})",
    colour = "Cases"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  transition_states(import_rate_label, transition_length = 2, state_length = 1) +
  ease_aes("cubic-in-out")

MMRcasemap <- animate(MMR_point_map, nframes = 100, fps = 10, width = 800, height = 600)
anim_save("MMRcumulative_cases.gif", MMRcasemap)
