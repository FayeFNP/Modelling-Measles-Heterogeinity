### DATA INTEGRATION AND PREPROCESSING -----------------------------------------

#Maps data from https://geoportal.statistics.gov.uk/

# Install necessary packages
if (!require(tidyverse)) {install.packages("tidyverse")}
if (!require(sf)) {install.packages("sf")}
if (!require(rgeoda)) {install.packages("rgeoda")}
if (!require(spdep)) {install.packages("spdep")}
if (!require(car)) {install.packages("car")}
if (!require(sfheaders)) {install.packages("sfheaders")}
if (!require(knitr)) {install.packages("knitr")}

# Load libraries
library(tidyverse)     
library(sf)
library(rgeoda)
library(spdep)
library(car)
library(sfheaders)
library(knitr)
library(broom)

ew_lads = st_read("./Upper_Tier_Local_Authorities_December_2022_Boundaries_UK_BFE_2393094386485431663/UTLA_MCTY_DEC_2022_UK_BFE.shp")
ew_lads %>% filter(substring(UTLA22CD,1,1)=='W') %>% select(UTLA22CD,UTLA22NM,-geometry) %>% print(n=22)

# Merge Rutland and Leicestershire
newrow <- ew_lads %>% filter(UTLA22NM=='Rutland' | UTLA22NM=='Leicestershire')
ew_lads <- ew_lads %>% filter(UTLA22NM!='Rutland' | UTLA22NM!='Leicestershire')

newrow <- newrow %>% 
  filter(UTLA22NM=='Leicestershire') %>% 
  mutate(geometry = 
           st_union(st_geometry(newrow %>% filter(UTLA22NM=='Leicestershire')),
                    st_geometry(newrow %>% filter(UTLA22NM=='Rutland'))))

ew_lads <- ew_lads %>% bind_rows(newrow)

vacc_data <- as_tibble(read.csv('./OneDrive_1_05-04-2025/child-vaccination-stats-csvs-2023-24/childhood-vaccination-la-num-denom-2023-24.csv', header=TRUE))


pop2022 <- readxl::read_xlsx('./OneDrive_1_05-04-2025/mye22tablesew2023geogs.xlsx',sheet='MYE3',skip=7) %>% 
  filter(is.element(Geography,c("Unitary Authority","London Borough","County","Metropolitan County")))


# Add Cumberland and Westmoreland and Furness to get Cumbria (old district)
newrow <- pop2022 %>% filter(is.element(Name,c('Cumberland','Westmorland and Furness'))) %>% 
  summarise(Code='E10000006',Name='Cumbria',Geography=unique(Geography),
            `Estimated Population mid-2021`=sum(`Estimated Population mid-2021`),
            Births=sum(Births),
            Deaths=sum(Deaths),
            `Births minus Deaths`=sum(`Births minus Deaths`),
            `Internal Migration Inflow`=sum(`Internal Migration Inflow`),
            `Internal Migration Outflow`=sum(`Internal Migration Outflow`),
            `Internal Migration Net`=sum(`Internal Migration Net`),
            `International Migration Outflow`=sum(`International Migration Outflow`),
            `International Migration Inflow`=sum(`Internal Migration Inflow`),
            `International Migration Net`=sum(`International Migration Net`),
            Other=sum(Other),
            `Estimated Population mid-2022`=sum(`Estimated Population mid-2022`))

pop2022 <- pop2022 %>% filter(!is.element(Name,c('Cumberland','Westmorland and Furness'))) %>% bind_rows(newrow)

# Patch Code to match shapefile (Change in boundary?) for North Yorkshire and Somerset
pop2022 <- pop2022 %>% mutate(Code=replace(Code, Code=='E06000065','E10000023'))
pop2022 <- pop2022 %>% mutate(Code=replace(Code, Code=='E06000066','E10000027'))

# England and Wales Local Authority Districts with matched population estimates
ew_lads <- ew_lads %>% inner_join(pop2022,by=c('UTLA22CD'='Code'))
ggplot(ew_lads) + geom_sf(aes(fill=log10(`Estimated Population mid-2022`)))

# Filter out proportion vaccinated by 5 years, 2 shots of MMR

MMRcoverage <- vacc_data %>% filter(Indicator == 'MMR2_5y')
denom <- vacc_data %>% filter(Indicator == '5y_Eligible_Pop') %>% select(Org_Code,n=Value)

MMRcoverage <- MMRcoverage %>% left_join(denom) %>% mutate(coverage = 100*Value/n)

MMRcoverage = as_tibble(read.csv('./MMRcoverage.csv'))

ggplot(ew_lads %>% inner_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))) + geom_sf(aes(fill=coverage))

ggplot(ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))) + geom_sf()

# Check local authorities that dont match to shape file
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))

# Check local authorities that don't match from shape file to MMR coverage
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))

dictionary <- as_tibble(read.csv('./MMR_shapefile_dictionary.csv'))
View(dictionary)


### MAPPING MMRCOVERAGE DATA TO SHAPE FILE --------------------------------------

# For Tyne and Wear
twrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Gateshead','Newcastle upon Tyne','North Tyneside','South Tyneside','Sunderland')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Gateshead','Newcastle upon Tyne','North Tyneside','South Tyneside','Sunderland')))
# Make new row
twrows <- twrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom)) %>% 
  mutate(`coverage`=100*Value/Denom)
twrows <- twrows %>% mutate(Org_Code='E11000007',Org_Name="Tyne and Wear",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(twrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For Greater Manchester
gmrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Bolton','Bury','Manchester','Oldham','Rochdale','Salford','Stockport','Tameside','Trafford','Wigan')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Bolton','Bury','Manchester','Oldham','Rochdale','Salford','Stockport','Tameside','Trafford','Wigan')))
# Make new row
gmrows <- gmrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom), .groups = 'drop') %>% 
  mutate(`coverage`=100*Value/Denom)
gmrows <- gmrows %>% mutate(Org_Code='E11000001',Org_Name="Greater Manchester",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(gmrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For Merseyside
mrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Knowsley','Liverpool','Sefton','St. Helens','Wirral')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Knowsley','Liverpool','Sefton','St. Helens','Wirral')))
# Make new row
mrows <- mrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom), .groups = 'drop') %>% 
  mutate(`coverage`=100*Value/Denom)
mrows <- mrows %>% mutate(Org_Code='E11000002',Org_Name="Merseyside",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(mrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For South Yorkshire
syrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Barnsley','Doncaster','Rotherham','Sheffield')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Barnsley','Doncaster','Rotherham','Sheffield')))
# Make new row
syrows <- syrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom), .groups = 'drop') %>% 
  mutate(`coverage`=100*Value/Denom)
syrows <- syrows %>% mutate(Org_Code='E11000003',Org_Name="South Yorkshire",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(syrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For West Yorkshire
wyrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Bradford','Calderdale','Kirklees','Leeds','Wakefield')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Bradford','Calderdale','Kirklees','Leeds','Wakefield')))
# Make new row
wyrows <- wyrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom), .groups = 'drop') %>% 
  mutate(`coverage`=100*Value/Denom)
wyrows <- wyrows %>% mutate(Org_Code='E11000006',Org_Name="West Yorkshire",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(wyrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For North Yorkshire - no aggregation needed
nyrow <- MMRcoverage %>% filter(Org_Name == 'North Yorkshire')
MMRcoverage <- MMRcoverage %>% filter(Org_Name != 'North Yorkshire')
# Update Org_Code
nyrow <- nyrow %>% mutate(Org_Code='E10000023')
# Add back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(nyrow)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For Somerset - no aggregation needed
srow <- MMRcoverage %>% filter(Org_Name == 'Somerset')
MMRcoverage <- MMRcoverage %>% filter(Org_Name != 'Somerset')
# Update Org_Code
srow <- srow %>% mutate(Org_Code='E10000027')
# Add back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(srow)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))


# For West Midlands
wmrows <- MMRcoverage %>% filter(is.element(Org_Name,c('Birmingham','Coventry','Dudley','Sandwell','Solihull','Walsall','Wolverhampton')))
MMRcoverage <- MMRcoverage %>% filter(!is.element(Org_Name,c('Birmingham','Coventry','Dudley','Sandwell','Solihull','Walsall','Wolverhampton')))
# Make new row
wmrows <- wmrows %>% 
  select(-`coverage`,-Org_Code,-Org_Name) %>% 
  group_by(CollectionYearRange,Parent_Org_Code,Parent_Org_Name,Child_Age,Indicator) %>% 
  summarise(Value=sum(Value),Denom=sum(Denom), .groups = 'drop') %>% 
  mutate(`coverage`=100*Value/Denom)
wmrows <- wmrows %>% mutate(Org_Code='E11000005',Org_Name="West Midlands",.before=Value)
# Add new aggregated row back into MMR Coverage table
MMRcoverage <- MMRcoverage %>% bind_rows(wmrows)
# Check matching
MMRcoverage %>% anti_join(ew_lads,by=c('Org_Code'='UTLA22CD'))
ew_lads %>% anti_join(MMRcoverage,by=c('UTLA22CD'='Org_Code'))

### MAPPING VACCINE COVERAGE----------------------------------------------------

# Join MMR coverage to shapefile data
MMR_map <- ew_lads %>% 
  inner_join(MMRcoverage, by = c("UTLA22CD" = "Org_Code"))

# Quick numeric summary
summary(MMR_map$coverage)

# Categorise into public health targets
MMR_map <- MMR_map %>% 
  mutate(threshold_group = factor(
    case_when(
      coverage >= 95 ~ "95%+ (Target Met)",
      coverage >= 90 ~ "90-94%",
      TRUE           ~ "<90% (Below Target)"
    ),
    levels = c("95%+ (Target Met)", "90-94%", "<90% (Below Target)")  
  ))

# Plot the map
ggplot(MMR_map) +
  geom_sf(aes(fill = threshold_group), colour = "grey30", size = .15) +
  scale_fill_manual(
    name   = "MMR Coverage Category",
    values = c(
      "95%+ (Target Met)"     = "#5F9E6E",  
      "90-94%"                = "#D89C30",  
      "<90% (Below Target)"   = "#8F4B57"   
    )
  ) +
  labs(
    title    = "MMR Vaccine Coverage by Local Authority (England & Wales)",
    subtitle = "Percentage of children with 2 MMR doses by age 5",
    caption  = "Data Source: UKHSA, ONS; https://geoportal.statistics.gov.uk/"
  ) +
  theme_minimal() 

# Breakdown of how many regions fall below certain thresholds
MMR_map %>%
  mutate(threshold_group = case_when(
    coverage >= 95 ~ "95%+ (Target Met)",
    coverage >= 90 ~ "90-94%",
    TRUE ~ "<90% (Below Target)"
  )) %>%
  count(threshold_group)

### BINOMIAL REGRESSION --------------------------------------------------------

names(MMR_map)

# Calculate predictors
MMR_map <- MMR_map %>%
  mutate(
    pop_2022 = `Estimated Population mid-2022`,
    area_km2 = as.numeric(st_area(geometry)) / 1e6,  # Area in square km
    pop_density = pop_2022 / area_km2,               # Population per km²
    birth_rate = Births / pop_2022,                  # Crude birth rate
    region = as.factor(Parent_Org_Name),             # Region as categorical variable
    int_mig_rate = 1000 * `Internal Migration Net` / pop_2022,        # per 1,000 people
    intl_mig_rate = 1000 * `International Migration Net` / pop_2022   # per 1,000 people
  )

# Fit the binomial regression model
model <- glm(
  formula = cbind(Value, Denom - Value) ~ log(pop_2022) + birth_rate + pop_density + int_mig_rate + intl_mig_rate + region,
  data = MMR_map,
  family = binomial(link = "logit")
)

# View model summary
summary(model)

# Convert to odds ratios
model_summary <- tidy(model, conf.int = TRUE, exponentiate = TRUE)

knitr::kable(model_summary, digits = 3, caption = "Binomial Regression Model for MMR Coverage")

#Predictive checks 

vif(model) # Check for Collinearity
step_model <- step(model, direction = "both") # Stepwise Regression
summary(step_model)

#Get predicted probabilities 
MMR_map$predicted <- predict(step_model, type = "response")
ggplot(MMR_map, aes(x = predicted, y = coverage / 100)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Predicted vs Observed MMR Coverage",
    x = "Predicted Coverage (Model)",
    y = "Observed Coverage (%)"
  ) +
  theme_minimal()

## Calculate McFadden’s Pseudo R²
1 - (step_model$deviance / step_model$null.deviance)

# plotting residuals
MMR_map$residuals <- residuals(step_model, type = "deviance")
ggplot(MMR_map, aes(x = predicted, y = residuals)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Residuals vs Predicted Coverage",
    x = "Predicted Coverage",
    y = "Deviance Residuals"
  ) +
  theme_minimal()

### HIERARCHICAL SPATIAL CLUSTERING --------------------------------------------

load('MMR_map.RData')

# Fill hole in Leicestershire that is breaking neighbour detection
MMR_map[142,] = sfheaders::sf_remove_holes(MMR_map[142,])
# Remove islands (Isle of Wight and Isle of Anglesey)
MMR_map <- MMR_map[-c(43,120),]

# Testing different number of clusters (2 to 20 clusters)
# For 2 clusters
library(sf)
w <- queen_weights(MMR_map) # Queen weights takes sf object not geoda
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(w, k = 2, MMR_map["coverage"])
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()

## For 3 clusters
w <- queen_weights(MMR_map) 
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(k=3,w=w,df = MMR_map %>% select(coverage))
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()

## For 6 clusters
w <- queen_weights(MMR_map) 
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(k=6,w=w,df = MMR_map %>% select(coverage))
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()

## For 10 clusters
w <- queen_weights(MMR_map)
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(k=10,w=w,df = MMR_map %>% select(coverage))
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)
checks <- spatial_validation(MMR_map,clusters$Clusters,w)
print(checks)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()

## For 15 clusters
w <- queen_weights(MMR_map)
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(k=15,w=w,df = MMR_map %>% select(coverage))
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()

## For 20 clusters
w <- queen_weights(MMR_map) 
any(!sf::st_is_valid(MMR_map))
sum(is.na(sf::st_geometry(MMR_map)))
clusters <- skater(k=20,w=w,df = MMR_map %>% select(coverage))
MMR_map$cluster_coverage <- as.factor(clusters$Clusters)

ggplot(MMR_map) +
  geom_sf(aes(fill = cluster_coverage)) +
  scale_fill_viridis_d(name = "Coverage Clusters") +
  labs(
    title = "Hierarchical Clusters of MMR Coverage",
    caption = "SKATER method on vaccine coverage only"
  ) +
  coord_sf(datum = NA) +
  theme_void()


#Clustering with multiple values of k and store RBTSS
k_range <- 2:20
wcss_values <- map_dbl(k_range, function(k) {
  print(k)
  result <- redcap(k = k, w = w, df = MMR_map %>% select(coverage))
  result$`The ratio of between to total sum of squares`  
})

#Plotting elbow curve
plot(k_range, wcss_values, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "The ratio of between to total sum of squares",
     main = "Elbow Plot for SKATER Clustering")


## Clustering with multiple values of k and store Entropy1
k_range <- 2:20
entropy_values <- map_dbl(k_range, function(k) {
  result <- skater(k = k, w = w, df = MMR_map %>% select(coverage))
  checks <- spatial_validation(MMR_map, result$Clusters, w)
  checks$Fragmentation$`Entropy*`  # Normalised entropy
})

plot(k_range, entropy_values, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Normalized Entropy",
     main = "Elbow Plot Using Entropy Measure")


## Clustering with multiple values of k and store Entropy2
k_range <- 2:20
entropy_values <- map_dbl(k_range, function(k) {
  result <- skater(k = k, w = w, df = MMR_map %>% select(coverage))
  checks <- spatial_validation(MMR_map, result$Clusters, w)
  checks$Fragmentation$Entropy # Normalised entropy
})

plot(k_range, entropy_values, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Normalized Entropy",
     main = "Elbow Plot Using Entropy Measure")

## Creating map of average MMR coverage per cluster
class(MMR_map)
dim(MMR_map)
head(MMR_map)

#Categorize coverage into defined % bands
MMR_map <- MMR_map %>%
  mutate(
    coverage_band = cut(
      coverage,
      breaks = c(60, 70, 80, 90, 95),
      include.lowest = TRUE,
      right = FALSE,
      labels = c("60–69%", "70–79%", "80–89%", "90–100%")
    )
  )

ggplot(MMR_map) +
  geom_sf(aes(fill = coverage_band)) +
  scale_fill_viridis_d(
    name = "MMR Coverage Band (%)",
    option = "C",
    direction = -1  # optional: reverse so high coverage is darker
  ) +
  labs(
    title = "Local Authorities Categorized by MMR Vaccine Coverage",
    subtitle = "Grouped in 10% intervals (age 5, 2 doses)",
    caption = "Data: UKHSA & ONS | Binned using fixed % thresholds"
  ) +
  coord_sf(datum = NA) +
  theme_void()

##MMR_map_backup <- MMR_map 

### RESIDUAL SPATIAL CLUSTERING ------------------------------------------------

#Global Morans I

nb <- poly2nb(MMR_map)
lw <- nb2listw(nb, style = "W")
global_moran_result <- moran.test(MMR_map$residuals, lw)
global_moran_result

#Local Morans I
local_moran <- localmoran(MMR_map$residuals, lw)

# Add the local Moran's I statistic to map
MMR_map$local_moran_i <- local_moran[,1]  # First column is I value

# Plot Map
ggplot(MMR_map) +
  geom_sf(aes(fill = local_moran_i)) +
  scale_fill_viridis_c(name = "Local Moran's I") +
  theme_minimal() +
  labs(title = "Local Spatial Clustering of Residuals (Local Moran's I)")

# Plot the improved map
ggplot(MMR_map) +
  geom_sf(aes(fill = local_moran_i)) +
  scale_fill_gradient2(
    low = "yellow",
    mid = "steelblue",
    high = "red",
    midpoint = 0,
    name = "Local Moran's I"
  ) +
  theme_void() +
  labs(
    title = "Local Spatial Clustering of Model Residuals",
    subtitle = "Red = High clustering, Blue = Low clustering",
    caption = "Local Moran's I statistic"
  )


## Plotting the LISA
colnames(local_moran)
local_moran <- localmoran(MMR_map$residuals, lw)
MMR_map$local_I <- local_moran[, "Ii"]
MMR_map$local_p <- local_moran[, "Pr(z != E(Ii))"]

# Standardize residuals (mean 0, sd 1)
MMR_map$std_resid <- scale(MMR_map$residuals)[, 1]

# Spatial lag of standardized residuals
lagged_resid <- lag.listw(lw, MMR_map$std_resid)

# Add lag to data
MMR_map$lagged_resid <- lagged_resid

# Classify cluster type
MMR_map$LISA_cluster <- "Not Significant"
MMR_map$LISA_cluster[MMR_map$std_resid > 0 & lagged_resid > 0 & MMR_map$local_p <= 0.05] <- "High-High"
MMR_map$LISA_cluster[MMR_map$std_resid < 0 & lagged_resid < 0 & MMR_map$local_p <= 0.05] <- "Low-Low"
MMR_map$LISA_cluster[MMR_map$std_resid > 0 & lagged_resid < 0 & MMR_map$local_p <= 0.05] <- "High-Low"
MMR_map$LISA_cluster[MMR_map$std_resid < 0 & lagged_resid > 0 & MMR_map$local_p <= 0.05] <- "Low-High"

# Convert to factor for plotting
MMR_map$LISA_cluster <- factor(MMR_map$LISA_cluster,
                               levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"))

ggplot(MMR_map) +
  geom_sf(aes(fill = LISA_cluster), color = "white", size = 0.1) +
  scale_fill_manual(
    name = "LISA Cluster Type",
    values = c(
      "High-High" = "#d7191c",   # Red
      "Low-Low" = "green",     # Green
      "High-Low" = "orange",    # Orange
      "Low-High" = "blue",    # Light Blue
      "Not Significant" = "steelblue"
    )
  ) +
  labs(
    title = "Local Moran's I Cluster Map",
    subtitle = "Spatial clustering of residuals from binomial regression",
    caption = "Only p ≤ 0.05 clusters shown as significant"
  ) +
  theme_void()

## Identifying areas of concern
MMR_map %>%
  filter(LISA_cluster %in% c("High-High", "Low-Low", "High-Low")) %>%
  select(UTLA22NM, LISA_cluster, coverage, residuals)

## Generating a table
MMR_map %>%
  st_drop_geometry() %>%  # Drop spatial geometry for a cleaner table
  filter(LISA_cluster %in% c("High-High", "Low-Low", "High-Low")) %>%
  select(`Local Authority` = UTLA22NM,
         `LISA Cluster` = LISA_cluster,
         `MMR Coverage (%)` = coverage,
         Residual = residuals) %>%
  arrange(`LISA Cluster`, desc(`MMR Coverage (%)`)) %>%
  kable(caption = "Local Authorities with Significant Local Spatial Clustering (LISA)")

### CITIES FOR SEEDING SIMULATION ----------------------------------------------

# Retrieve regression model
MMR_map$predicted <- predict(step_model, type = "response")

# Identifying areas with lowest predicted MMR coverage:
MMR_map %>%
  st_drop_geometry() %>%
  select(UTLA22NM, predicted) %>%
  mutate(predicted = round(predicted * 100, 2)) %>%
  arrange(predicted) %>%
  head(50)

low_coverage <- MMR_map %>%
  st_drop_geometry() %>%
  select(UTLA22NM, predicted) %>%
  mutate(predicted = round(predicted * 100, 2)) %>%
  arrange(predicted) %>%
  head(50)

write.csv(low_coverage, "low_coverage_districts.csv", row.names = FALSE)

# Identifying areas with highest predicted MMR coverage:
MMR_map %>%
  st_drop_geometry() %>%
  select(UTLA22NM, predicted) %>%
  mutate(predicted = round(predicted * 100, 2)) %>%
  arrange(desc(predicted)) %>%
  head(10)
